## Directed evolution simulation analysis
# August/20 v2.1

#= While trying different measures, I found that the KL divergence 
of the conditional probabilities - comparing r1 to r0 gave the most intuitive results.
As such, I implement it here for further testing =#

using BenchmarkTools
using Plots
using SharedArrays
using DelimitedFiles
using StatsBase

## Dependencies from selection file:
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H' , 'I', 'K', 'L', 'M']#, 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # Not using a 'packaged' amino acid symbol to allow unnatural amino acids and deletions to be included

## Interesting variables to pass between files:
    # fitness_function

## Importing data from selection
    r0_pop = readdlm("r0_pop.csv", ',', Char)
    r0_pop = convert(SharedArray, r0_pop)

    r1_pop = readdlm("r1_pop.csv", ',', Char)
    r1_pop = convert(SharedArray, r1_pop)

    sequence_length = size(r0_pop, 2)


## Functions
    function z_subsets(population::SharedArray{Char, 2},
            amino_acids::Array{Char, 1}, 
            sequence_length::Int64)

        # This returns an array containing the indexes for each
        # amino_acid in a sorted array (x for residues, y for aa)

        indexes = zeros(Int64, length(amino_acids), sequence_length)

        for z = 1:sequence_length
            indexes[:, z] = z_subsets_aa(population, z, amino_acids, sequence_length)
        end
        return indexes
    end

    function z_subsets_aa(population::SharedArray{Char, 2},
                z::Int64, amino_acids::Array{Char, 1}, 
                sequence_length::Int64)

        # This returns an array containing the indexes for each
        # amino_acid in a sorted array (x for residues, y for aa) for a given amino acid (z)

        indexes = zeros(Int64, length(amino_acids))
        pop_sort = sortslices(population, dims = 1, by= a -> a[z])

        for aa = 1: length(amino_acids)
            if aa == 1 && findlast(b->b == amino_acids[aa], pop_sort[:,z]) === nothing
                indexes[aa] = 1
            elseif findlast(b->b == amino_acids[aa], pop_sort[:,z]) === nothing
                indexes[aa] = indexes[aa-1]
            else
                indexes[aa] = findlast(b->b == amino_acids[aa], pop_sort[:,z]) 
            end
        end
        return indexes
    end

    function givenz_subset_freqs(population::SharedArray{Char, 2},
                residue::Int64,
                amino_acids::Array{Char, 1})

        sequence_length = size(population, 2)
        population_index = z_subsets(population, amino_acids, sequence_length)

        population_sorted = sortslices(population, dims = 1, by= res -> res[residue])
        subpop_freqs = zeros(Float64, length(amino_acids), sequence_length , length(amino_acids))

        for aa = 1:length(amino_acids)
            if aa == 1
                first_index = 1
                final_index = population_index[aa, residue]
                subset = convert(SharedArray, population_sorted[first_index:final_index, :])
                subpop = z_subsets(subset, amino_acids, sequence_length)
            else
                first_index = population_index[aa-1, residue] + 1
                final_index = population_index[aa, residue]
                subset = convert(SharedArray, population_sorted[first_index:final_index, :])
                subpop = z_subsets(subset, amino_acids, sequence_length)
            end
        subpop_freqs[:,:,aa] = z_frequency(subpop)
        end

        return subpop_freqs
    end

    function z_frequency(subset_index::Array{Int64, 2})
        # this function converts the indexes identified in z_subsets into a per residue frequency table 
        subset_freq = zeros(Float64, size(subset_index, 1), size(subset_index, 2))
    
        for res = 1 : sequence_length
            for a = 1: size(subset_index, 1)
                if a == 1
                    freq = (subset_index[a, res] - 1)/(subset_index[size(subset_index, 1), res])
                    subset_freq[a, res] = freq 
                else
                    freq = (subset_index[a, res] - subset_index[a-1, res])/(subset_index[size(subset_index, 1), res])
                    subset_freq[a, res] = freq 
                end
            end
        end
    
        return subset_freq
    end

    function kulleib_1d_aa(r1_freq::Array{Float64, 2},
                r0_freq::Array{Float64, 2},
                residue::Int64)

        # this function returns the Kullback_Leibler Divergence
        # for an individual residue

        kl = 0.0
        for a = 1:size(r1_freq, 1)
            if r0_freq[a, residue] != 0

            # ISSUE: where r0_freq=0, log2 -> Inf. These are being ignored from the calculation.
            # Not sure if that is mathematically sound

                res_contrib = r1_freq[a, residue] * (log2(r1_freq[a, residue]) - log2(r0_freq[a, residue]))
                kl += res_contrib
            end
        end
        return kl
        # Positive values mean sequence convergence (i.e. selection)
        # Negative values mean sequence divergence
    end

    function kulleib_1d(r1_freq::Array{Float64, 2},
            r0_freq::Array{Float64, 2})

        # this function returns the Kullback_Leibler Divergence (r1||r0)
        # for all individual residues in a sequence. This function can be
        # extended to compare against an environment-specific r0 if 
        # phylogenetic data is being used as r1

        kl = zeros(Float64, size(r1_freq, 2), 1)
        for residue = 1: size(r1_freq, 2)
            kl[residue] = kulleib_1d_aa(r1_freq, r0_freq, residue)
        end

        return kl

    end



## Mapping subsets
    r0_index = z_subsets(r0_pop, amino_acids, sequence_length)
    r1_index = z_subsets(r1_pop, amino_acids, sequence_length)

    ## Population given a fixed residue
    given_res_r0 = Array{Array{Float64,3},1}(undef,sequence_length)
    given_res_r1 = Array{Array{Float64,3},1}(undef,sequence_length)

    for res = 1: sequence_length
        given_res_r0[res] = givenz_subset_freqs(r0_pop, res, amino_acids)
        given_res_r1[res] = givenz_subset_freqs(r1_pop, res, amino_acids)
    end 
    
    ## KL given a fixed residue
    given_kla = zeros(Float64, sequence_length, length(amino_acids), sequence_length)

    for given = 1: sequence_length
        for aa = 1: length(amino_acids)
            given_kla[:,aa,given] = kulleib_1d(given_res_r1[given][:,:,aa], given_res_r0[given][:,:,aa])
        end
    end
    
    


## Showing results for all 
    plot1 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 1", ylabel = "Other_residues", given_kla[:,:,1]);
    plot2 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 2", ylabel = "Other_residues", given_kla[:,:,2]);
    plot3 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 3", ylabel = "Other_residues", given_kla[:,:,3]);
    plot4 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 4", ylabel = "Other_residues", given_kla[:,:,4]);
    plot5 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 5", ylabel = "Other_residues", given_kla[:,:,5]);
    plot6 = heatmap([string(aa) for aa in amino_acids],[string(res) for res = 1:sequence_length],
        xlabel = "Residue 6", ylabel = "Other_residues", given_kla[:,:,6]);

    plot(plot1, plot2, plot3, plot4, plot5, plot6,  clims=(0,0.2), layout = (3,2))

## Histogram
    simple_kl_givenz = reshape(given_kla, (length(given_kla)))
    simple_kl_givenz = sort(simple_kl_givenz)

    plot([simple_kl_givenz], bins = 20, seriestype = [:barhist],
    fillcolor =[:lightblue], linecolor=[:black],thickness_scaling = 1,
    label = "KL divergence (r1|z || r0|z)")     

    println("The following interactions may be of importance:")
    for hit = 1:10
        combinations = findlast(x->x==simple_kl_givenz[end+1-hit], given_kla)
        println("$hit: Residue $(combinations[3]) when $(amino_acids[combinations[2]]) is linked to residue $(combinations[1])-> $(simple_kl_givenz[end+1-hit])") 
    end

## Inspecting the distributions per residue, given a fixed residue
fixed = 3       # residue number
fixed_aa = 3    # amino acid # in amino_acids
res = 2         # residue being looked at

sel_effect_givenz = plot(bar([string(aa) for aa in amino_acids], given_res_r0[fixed][:, res, fixed_aa]),
bar([string(aa) for aa in amino_acids],given_res_r1[fixed][:, res, fixed_aa]),
fillcolor =[:lightgrey :blue], ylims =(0,1), label = ["R0_givenz" "R1_givenz"])
## Directed evolution simulation analysis
# July/20 v3.1


# The idea is to create a simulation of directed evolution in which:
# - given a set of parameters, a new population can be created (R0)
# - given a set of stochastic criteria, R0 can be selected
# - The resulting population (R1) includes a component of noise
# - The program outputs key information of the run.

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


## functions - 1D analysis
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
    z::Int64,
    amino_acids::Array{Char, 1}, 
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

## functions - 2D analysis
# MI original implementation
function probabilities(x)
    counts = countmap(x)
    probs = values(counts)./sum(values(counts))
    return probs
end
    
function mi_xy2d(r_pop::SharedArray{Char, 2})
    
    ro_r, co_r = size(r_pop)
    information_xixj_r = zeros((co_r, co_r))
    
    for i in 1:co_r
        xi= r_pop[:,i]
        for j in range(1, length=co_r)
            xj= r_pop[:,j]
            xixj   = collect(zip(xi, xj))
            pxi      = probabilities(xi)
            pxj      = probabilities(xj)
            pxixj    = probabilities(xixj)
            entropy_xi     = (-1.0).*sum(log2.(pxi).*pxi)
            entropy_xj     = (-1.0).*sum(log2.(pxj).*pxj)
            entropy_xixj   = (-1.0).*sum(log2.(pxixj).*pxixj)
            information_xixj_r[i,j] = entropy_xi + entropy_xj - entropy_xixj
        end
    end
    
    return information_xixj_r
end

# MI fresh implementation
# MI = SUMx SUMy p(x, y) * log[p(x,y)/(p(x)*p(y))]

function mi_2points(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1},
    x::Int64, y::Int64)

    prob_pair = zeros(Float64, length(amino_acids), length(amino_acids))

    for call = 1: size(population, 1)
        res1 = population[call, :][x]
        res2 = population[call, :][y] 
        prob_pair[findfirst(aa->aa==res1, amino_acids)[1], findfirst(aa->aa==res2, amino_acids)[1]] += 1.0
    end

    prob_pair = prob_pair/sum(prob_pair)

    mi = 0

    for a = 1:length(amino_acids)
        py = sum(prob_pair[a,:])
        px = sum(prob_pair[:,a])
        pxy = prob_pair[a,a]
        if pxy != 0
            mi += pxy * log2((pxy)/(px * py))
        end
    end

    return mi
end

function mi_xy2d2(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1})

    mi_2d = zeros(Float64, size(population, 2), size(population, 2))
    for res1 = 1: size(population, 2)
        for res2 = 1: size(population, 2)
            mi_2d[res1,res2] = mi_2points(population, amino_acids, res2, res1)
        end
    end
    return mi_2d
end

# Kullback-Leibler Divergence on the conditional probability - subset compared to population

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

function kulleib_1d_givenz(given_res_freq::Array{Float64, 3}, 
    population_freq::Array{Float64, 2})

    kl = zeros(Float64, size(population_freq, 1), size(population_freq, 2))
    for aa = 1: size(given_res_freq, 3)
        kl[aa,:] = kulleib_1d(given_res_freq[:,:,aa], population_freq)
    end
    return kl
end

## functions - 3D analysis
# Conditional probability
# p(x,y|z)

function freq_xy(population::SharedArray{Char, 2},
    x::Int64,
    y::Int64,
    amino_acids::Array{Char, 1})

    # Creates a frequency table comparing positions x and y
    freq_map = zeros(length(amino_acids), length(amino_acids))

    if size(population, 1) == 0
        return freq_map
    else
        for entry = 1: size(population, 1)
            xi = findall(a->a== population[entry,x], amino_acids)
            yi = findall(b->b== population[entry,y], amino_acids)
            freq_map[yi[1], xi[1]] += 1
        end
    end

    freq_map = freq_map ./ sum(freq_map)
    return freq_map
end

function pxy_givenz(population::SharedArray{Char, 2}, 
    xy::Array{Int64,1}, z::Int64,  
    amino_acids::Array{Char, 1})

    pop_sort = sortslices(population, dims = 1, by= a -> a[z])
    indexes = z_subsets_aa(population, z, amino_acids, sequence_length)

    pxiyi_givenz = zeros(Float64, length(amino_acids), length(amino_acids), length(amino_acids))
    for zi = 1:length(amino_acids)
        if zi == 1
            subset = convert(SharedArray{Char, 2},pop_sort[1:indexes[zi], :])
            pxiyi_givenz[:,:,zi] = freq_xy(subset, xy[1], xy[2], amino_acids)
        else zi == 1
            subset = convert(SharedArray{Char, 2},pop_sort[(indexes[zi-1]+1):indexes[zi], :])
            pxiyi_givenz[:,:,zi] = freq_xy(subset, xy[1], xy[2], amino_acids)
        end
    end

    return pxiyi_givenz
end

function delta_pxy_givenz(r1_population::SharedArray{Char, 2}, 
    r0_population::SharedArray{Char, 2}, 
    xy::Array{Int64,1}, z::Int64, 
    amino_acids::Array{Char, 1})

    delta_p = zeros(Float64, length(amino_acids),length(amino_acids),length(amino_acids))
    delta_p = pxy_givenz(r1_population, xy ,z, amino_acids) - pxy_givenz(r0_population, xy,z, amino_acids)

    return delta_p
end

# MI original implementation
# I(X;Y;Z)

function mi_xyz3d(r_pop::SharedArray{Char, 2})

    ro_r, co_r = size(r_pop);
    information_xixjxk_r = zeros((co_r, co_r, co_r));
    for i in 1:co_r;
        xi= r_pop[:,i];
        for j in range(1, length=co_r);
            xj= r_pop[:,j];
            for k in range(1, length=co_r);
                xk= r_pop[:,k];
                xixj   = collect(zip(xi, xj));
                xixk   = collect(zip(xi, xk));
                xjxk   = collect(zip(xj, xk));
                xixjxk = collect(zip(xi, xj, xk));

                pxi      = probabilities(xi);
                pxj      = probabilities(xj);
                pxk      = probabilities(xk);
                pxixj    = probabilities(xixj);
                pxixk    = probabilities(xixk);
                pxjxk    = probabilities(xjxk);
                pxixjxk  = probabilities(xixjxk);

                entropy_xi     = (-1.0).*sum(log2.(pxi).*pxi);
                entropy_xj     = (-1.0).*sum(log2.(pxj).*pxj);
                entropy_xk     = (-1.0).*sum(log2.(pxk).*pxk);
                entropy_xixj   = (-1.0).*sum(log2.(pxixj).*pxixj);
                entropy_xixk   = (-1.0).*sum(log2.(pxixk).*pxixk);
                entropy_xjxk   = (-1.0).*sum(log2.(pxjxk).*pxjxk);
                entropy_xixjxk = (-1.0).*sum(log2.(pxixjxk).*pxixjxk);

                information_xixjxk_r[i,j,k]    = (entropy_xi + entropy_xj + entropy_xk) -
                    (entropy_xixj + entropy_xixk + entropy_xjxk) + (entropy_xixjxk);
            end
        end
    end

    return information_xixjxk_r
end

# Interaction information or conditional mutual information - I(X;Y|Z)
# SUMx SUMy SUMz p(x, y, z) * log[(p(z)*p(x,y,z))/(p(x,z)*p(y,z))]

function info_xy_z(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1},
    x::Int64, y::Int64, z::Int64)

    prob_trio = zeros(Float64, length(amino_acids), length(amino_acids), length(amino_acids))

    for call = 1: size(population, 1)
        res1 = population[call, :][x]
        res2 = population[call, :][y]
        res3 = population[call, :][z]
        prob_trio[findall(aa->aa==res1, amino_acids)[1], findall(aa->aa==res2, amino_acids)[1],
        findall(aa->aa==res3, amino_acids)[1]] += 1.0
    end

    prob_trio = prob_trio/sum(prob_trio)

    info = 0

    for a = 1:length(amino_acids)
        pz   = sum(prob_trio[:,:,a])
        pxyz = prob_trio[a,a,a]
        pxz  = sum(prob_trio[a,:,a]) 
        pyz  = sum(prob_trio[:,a,a])
        if pxyz != 0
            info += pxyz * log2((pz*pxyz)/(pxz*pyz))
        
        end
    end

    return info
end

function info_xy_givenz(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1})

    mi_3d = zeros(Float64, size(population, 2), size(population, 2),size(population, 2))
    for res1 = 1: size(population, 2)
        for res2 = 1: size(population, 2)
            for res3 = 1: size(population, 2)
                mi_3d[res1,res2,res3] = info_xy_z(population, amino_acids, res2, res1, res3)
            end
        end
    end
    return mi_3d
end

# Specific Interaction information - De Cruys SI1
# log[(p(x,y)*p(y,z)*p(x,z))/(p(x)*p(y)*p(z)*p(xyz)] - original file I shared had the wrong calculation.
# However, I fixed the calculation and the result was meaningless.

## Total correlation - De Cruys SI2
# log[(p(x,y,z)/(p(x)*p(y)*p(z)] - calculation is incorrect


## 1D analysis - Kullback_Leibler divergence
    ## Creating subsets
    r0_index = z_subsets(r0_pop, amino_acids, sequence_length)
    r1_index = z_subsets(r1_pop, amino_acids, sequence_length)
    ## Converting those subsets into frequencies
    r0_freq = z_frequency(r0_index)
    r1_freq = z_frequency(r1_index)
    ## Inspecting the distributions per residue to compare pre- and post-selection
    z = 1
    sel_effect = plot(bar(r0_freq[:, z]), bar(r1_freq[:, z]),
    fillcolor =[:lightgrey :blue], ylims =(0,1), label = ["R0" "R1"])
    ## Kullback-Leibler Divergence for individual residues
    single_kl = kulleib_1d(r1_freq, r0_freq)
    kl_seq = plot(bar(single_kl), fillcolor =[:lightgrey], label = "KL divergence (r1||r0)") 

## 2D analysis (original MI)
    ## Change in mutual information (heatmap)
    r0_mi2d = mi_xy2d(r0_pop)
    r1_mi2d = mi_xy2d(r1_pop)

    delta_mi2d = r0_mi2d - r1_mi2d


    plot(heatmap(r0_mi2d),heatmap(r1_mi2d), heatmap(delta_mi2d),
    title = ["R0_MI" "R1_MI" "delta_MI"])
    ## Change in mutual information (distribution)
    simple_deltami2d = reshape(delta_mi2d, (length(delta_mi2d)))

    plot([simple_deltami2d], bins = 20, seriestype = [:barhist],
    fillcolor =[:orange], linecolor=[:black],thickness_scaling = 1,
    label = "deltaMI 2D")
    
## 2D analysis (re-implemented MI)
    ## Change in mutual information (heatmap)
    r0_mi2d2 = mi_xy2d2(r0_pop, amino_acids)
    # getting an error in this function here but not on notebook
    r1_mi2d2 = mi_xy2d2(r1_pop, amino_acids)

    delta_mi2d2 = r0_mi2d2 - r1_mi2d2


    plot(heatmap(r0_mi2d2),heatmap(r1_mi2d2), heatmap(delta_mi2d2),
    title = ["R0_MI" "R1_MI" "delta_MI"])

    simple_deltami2d2 = reshape(delta_mi2d2, (length(delta_mi2d2)))

    plot([simple_deltami2d2], bins = 20, seriestype = [:barhist],
    fillcolor =[:purple], linecolor=[:black],thickness_scaling = 1,
    label = "deltaMI 2D") 

## Comparison between the 2 MI calculations
    plot([simple_deltami2d, simple_deltami2d2], bins = 20, seriestype = [:barhist :stephist],
    fillcolor =[:orange :purple], linecolor=[:orange :purple],thickness_scaling = 2,
    label = ["sum of H" "p(x)"]) 

#= Although the graphs are not superimposable, I think the differences may simply be
the result of approximation errors by either or both algorithms used =#

## 2D analysis (Kullback-Leibler Divergence on the conditional probability - subset compared to population)
    ## Mapping subsets
    r1_index = z_subsets(r1_pop, amino_acids, sequence_length)

    ## Average population
    r1_freq = z_frequency(r1_index)

    ## Population given a fixed residue
    given_res = Array{Array{Float64,3},1}(undef,sequence_length)
    for res = 1: sequence_length
        given_res[res] = givenz_subset_freqs(r1_pop, res, amino_acids)
        # this array of tensors contain all subpopulations
    end

    ## KL given a fixed residue
    given_kl = zeros(Float64, length(amino_acids), sequence_length, sequence_length)
    for res = 1: sequence_length
        given_kl[:,:,res] = kulleib_1d_givenz(given_res[res], r1_freq)
        # this array of tensors contain all possible KL divergences (r1|z || r1)
    end

    ## Inspecting the distributions per residue, given a fixed residue
    fixed = 4
    fixed_aa = 4
    res = 3
    sel_effect_givenz = plot(bar([string(aa) for aa in amino_acids], r1_freq[:, res]), bar([string(aa) for aa in amino_acids],
        given_res[fixed][:, res, fixed_aa]), fillcolor =[:lightgrey :blue], ylims =(0,1), label = ["R1" "R1_givenz"])

    res = 3
    plot(heatmap([string(res) for res = 1:sequence_length],[string(aa) for aa in amino_acids], ylabel = "Residue $res",
        xlabel = "Other_residues", title = "KL divergence (r1|z || r1)", given_kl[:,:,res]))

## 2D analysis (Kullback-Leibler Divergence on the conditional probability - r1 subset compared to r0 subset)

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
    
    ## Inspecting the distributions per residue, given a fixed residue
    fixed = 3
    fixed_aa = 4
    res = 4

    sel_effect_givenz = plot(bar([string(aa) for aa in amino_acids], given_res_r0[fixed][:, res, fixed_aa]),
    bar([string(aa) for aa in amino_acids],given_res_r1[fixed][:, res, fixed_aa]),
    fillcolor =[:lightgrey :blue], ylims =(0,1), label = ["R0_givenz" "R1_givenz"])

    res = 3
    plot(heatmap([string(aa) for aa in amino_acids],
    [string(res) for res = 1:sequence_length], xlabel = "Residue $res",
    ylabel = "Other_residues", title = "KL divergence (r1|z || r0|z)", clims=(0,0.2),
    given_kla[:,:,res]))

    simple_kl_givenz = reshape(given_kla, (length(given_kla)))

    plot([simple_kl_givenz], bins = 20, seriestype = [:barhist],
    fillcolor =[:lightblue], linecolor=[:black],thickness_scaling = 1,
    label = "KL divergence (r1|z || r0|z)") 
     
## Changes in conditional probability
    xy = [2,3]
    z = 4
    deltapxy_givenzaa = delta_pxy_givenz(r1_pop, r0_pop, xy, z, amino_acids);

    aa = 3
    plot(heatmap([string(aa) for aa in amino_acids],[string(aa) for aa in amino_acids], deltapxy_givenzaa[:,:,aa]), clims = (minimum(deltapxy_givenzaa), maximum(deltapxy_givenzaa)),
    xlabel = "Residue $(xy[1])", ylabel = "Residue $(xy[2])", title = "Given residue $z as $(amino_acids[aa])")

    ## Simple analysis of delta p(x,y)|z
    simple_deltapxy = reshape(deltapxy_givenzaa, (length(deltapxy_givenzaa)))
    simple_deltapxy = sort(simple_deltapxy)
    plot([simple_deltapxy], bins = simple_deltapxy[1]:0.001:simple_deltapxy[end], seriestype = [:barhist],
    fillcolor =[:lightgreen], linecolor=[:black],thickness_scaling = 1.5, label = "delta p(x,y|z)")
    
## 3D analysis (MI)
    ## Delta MI in 3D

    r0_mi3d = mi_xyz3d(r0_pop)
    r1_mi3d = mi_xyz3d(r1_pop)

    delta_mi3d = r0_mi3d - r1_mi3d;
    
    aa = 2
    plot(heatmap([string(res) for res = 1:sequence_length],
    [string(res) for res = 1:sequence_length], xlabel = "Residue",
    ylabel = "Residue", title = "deltaMI with residue $aa", clims=(0,0.2),
    delta_mi3d[:,:,aa]))

    ## Simple analysis of MI in 3D
    simple_deltami3d = reshape(delta_mi3d, (length(delta_mi3d)))
    simple_deltami3d = sort(simple_deltami3d);

    plot([simple_deltami3d], bins = simple_deltami3d[begin]:0.01:simple_deltami3d[end],
            seriestype = [:barhist], fillcolor =[:orange], linecolor=[:orange],thickness_scaling = 1,
            label = ["deltaMI 3D"]) 

## 3D analysis - Interaction information or conditional mutual information - I(X;Y|Z)
    # SUMx SUMy SUMz p(x, y, z) * log[(p(z)*p(x,y,z))/(p(x,z)*p(y,z))]   
    r0_infoxy = info_xy_givenz(r0_pop, amino_acids)
    r1_infoxy = info_xy_givenz(r1_pop, amino_acids)

    delta_infoxy = r0_infoxy - r1_infoxy;

    aa = 3
    plot(heatmap([string(res) for res = 1:sequence_length],
    [string(res) for res = 1:sequence_length], xlabel = "Residue",
    ylabel = "Residue", title = "I(X;Y|Z) for Z = Residue $aa", clims=(-0.05,0.05),c=cgrad([:blue, :white, :red]),
    delta_infoxy[:,:,aa]))
    
    # Simple analysis of MI in 3D
    simple_deltainfoxy = reshape(delta_infoxy, (length(delta_infoxy)))
    simple_deltainfoxy = sort(simple_deltainfoxy);
    plot([simple_deltainfoxy], bins = simple_deltainfoxy[begin]:0.002:simple_deltainfoxy[end],
            seriestype = [:barhist], fillcolor =[:orange], linecolor=[:orange],thickness_scaling = 1,
            label = ["I(X;Y|Z)"]) 

    
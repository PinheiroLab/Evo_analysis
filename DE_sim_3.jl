## Directed evolution simulation and analysis
# July/20 v3.0


# The idea is to create a simulation of directed evolution in which:
# - given a set of parameters, a new population can be created (R0)
# - given a set of stochastic criteria, R0 can be selected
# - The resulting population (R1) includes a component of noise

using BenchmarkTools
using Plots
#using BioSequences
#using FreqTables
#using DataFrames
using Plots
#using Bio.Seq
using StatsBase
using Distributed
using SharedArrays
#using DelimitedFiles


## Library generation functions
function position!(residue::Int64,
                   allowed_aa::Array{Char,1},
                   position_composition::Array{Pair{Int64,Array{Char,1}},1})
    # the function sets the sequence constraints for a given position in a synthesised library

    new_pos = Pair(residue, allowed_aa)
    return vcat(position_composition, new_pos)
end

function r0_library(position_composition::Array{Pair{Int64,Array{Char,1}},1},
                    r0_pop_size::Int64,
                    sequence_length::Int64)
    # the function generates a library of sequences based on the positional_composition requested

    r0_pop = SharedArray{Char, 2}(r0_pop_size, sequence_length)
    @distributed for position = 1:length(position_composition)
        for entry = 1:r0_pop_size
            residue = position_composition[position][1]
            pick = rand(position_composition[position][2])
            r0_pop[entry, residue] = pick
        end
    end
    return r0_pop
end

function pos_criteria_check(position_composition::Array{Pair{Int64,Array{Char,1}},1},
                            sequence_length::Int64,
                            amino_acids::Array{Char,1})
    # this function checks that all positions are accounted for in the design of the R0_library
    # non-declared positions become randomised

    max_conditions = 1:sequence_length
    for condition = 1:length(position_composition)
        if position_composition[condition][1] > sequence_length
            error("Residue distribution outside of the sequence length")
        elseif position_composition[condition][1] ⊈ max_conditions
            error("Residue distribution may have a duplication")
        else
            max_conditions = setdiff(max_conditions, position_composition[condition][1])
        end
    end

    if max_conditions != []
        for remain = 1:length(max_conditions)
            position_composition = position!(max_conditions[remain], amino_acids, position_composition)
        end
    end
    return sort!(position_composition, by = x -> x[1])
end

## Library fitness generation functions
function fit_criteria!(residues::Array{Int64, 1},
                       high_fit_combinations::Array{Array{Char,1},1},
                       fitness_function::Array{Pair{Array{Int64,1},Array{Array{Char,1},1}},1})
    # the function helps create a list of fitness criteria

    fit_comb = Pair(residues, high_fit_combinations)
    return vcat(fitness_function, fit_comb)
end

function fitness_seq(sequence::Array{Char, 1},
                     fitness_function::Array{Pair{Array{Int64,1},Array{Array{Char,1},1}},1},
                     fit_states::Array{Float64, 1})

    # this function sums up the fitness of a single sequence
    fitness = 0.0
    for criteria = 1:length(fitness_function)
        if sequence[fitness_function[criteria][1]] in fitness_function[criteria][2]
            fitness += fit_states[end]
        else
            fitness += fit_states[1]
            end
    end
    return fitness
end

function fitness_calculation(r0_pop::SharedArray{Char, 2},
        fitness_states::Array{Float64, 1},
        fitness_function::Array{Pair{Array{Int64,1},Array{Array{Char,1},1}},1})
    # this function calculates the individual sequence fitness for a population

    #= I'm creating race conditions in this function when I paralelize. I will
    need some thinking here to re-write the function to stop that. =#

    fitness = SharedArray{Float32, 1}(size(r0_pop, 1), init = 0.0)     # Creates an empty fitness matrix
    for entry = 1:size(r0_pop,1)
        sequences = r0_pop[entry, :]
        fitness[entry,1] = fitness_seq(sequences, fitness_function, fitness_states)
    end
    return fitness

end

## Selection and noise functions

function selection(r0_pop::SharedArray{Char, 2},
                   fitness::SharedArray{Float32,1},
                   trials::Int64,
                   sequence_length::Int64)

    r1_pop = SharedArray{Char, 2}(0, sequence_length)

    #= I'm creating race conditions in this function when I paralelize. I will
    need some thinking here to re-write the function to stop that. =#

    for trial = 1:trials
        pick = rand(1:size(r0_pop,1))
        p_recovery = fitness[pick]
        if p_recovery >= rand()
            pick_seq = reshape(r0_pop[pick,:],(1,sequence_length))
            r1_pop = vcat(pick_seq, r1_pop)
        end
    end
    return r1_pop
end

function sel_noise(noise_level::Float64,
                    r0_pop::SharedArray{Char,2},
                    sequence_length::Int64)

    #= I'm creating race conditions in this function when I paralelize. I will
    need some thinking here to re-write the function to stop that. =#

    noise_pop = Array{Char, 2}(undef, 0, sequence_length)
    noise_pop_size = Int64(round(noise_level*size(r0_pop, 1)))
    for noise = 1: noise_pop_size
        pick = rand(1:size(r0_pop,1))
        pick_seq = reshape(r0_pop[pick,:],(1,sequence_length))
        noise_pop = vcat(pick_seq, noise_pop)
    end
    return noise_pop
end

function noise_amp(population::SharedArray{Char,2},
                    end_size::Int64)

    #= I'm creating race conditions in this function when I paralelize. I will
    need some thinking here to re-write the function to stop that. =#

    r1_pop = population
    if size(population, 1) < end_size
        amplification_noise = Int64(end_size - size(population, 1))
        r1_pop = SharedArray{Char,2}(end_size, size(population, 2))

        for transfer = 1:size(population, 1)
            r1_pop[transfer, :] = population[transfer, :]
        end

        for amplif = 1: amplification_noise
            pick = rand(1:(size(population, 1)+amplif-1))
            pick_seq = reshape(r1_pop[pick,:],(1,size(population,2)))
            r1_pop[(size(population, 1)+amplif), :] = pick_seq
        end
    end
    return r1_pop
end

function mutation!(population::SharedArray{Char,2},
                    error_rate::Float64,
                    amino_acids::Array{Char, 1})

    @distributed for position = 1:length(population)
        if error_rate >= rand()
            population[position] = rand(amino_acids)
        end
    end
    return population
end

## General constraints
    r0_pop_size = Int64(1e5)     # R0 size
    sequence_length = 6          # Seqence Length
    trials = Int64(5e5)          # Number of selection events.
    r1_pop_size = Int64(1e5)     # Size of the end population

    #= The final population can exceed r1_pop_size if selected sequences and noise
    exceed this limit. Adjust parameters of the run so that this does not occur =#

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M']#, 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # Not using a 'packaged' amino acid symbol to allow unnatural amino acids and deletions to be included

    ## Creating R0 - the starting population
    position_composition = Pair{Int64,Array{Char,1}}[]
    position_composition = position!(1, ['A'], position_composition)
    #position_composition = position!(2, ['H'], position_composition)
    position_composition = pos_criteria_check(position_composition, sequence_length, amino_acids)
    # This sets the positional requirements in R0
    r0_pop = r0_library(position_composition, r0_pop_size, sequence_length)
    # This generates the random R0

## Fitness landscape functions
    fitness_function = Pair{Array{Int64,1},Array{Array{Char,1},1}}[]
    # Resets the functional landscape
    #fitness_function = fit_criteria!([5, 6, 7], [['A', 'A', 'A'],['D', 'D', 'D'],['D', 'G', 'G']], fitness_function)
    fitness_function = fit_criteria!([1, 3], [['A', 'C'],['A', 'D'],['H', 'H']], fitness_function)
    fitness_function = fit_criteria!([2], [['H'],['A']], fitness_function)
    #fitness_function = fit_criteria!([3], [['A']], fitness_function)
    #fitness_function = fit_criteria!([4], [['A']], fitness_function)

    # This creates a fitness function
    fit_lo = 0.0005                     # probability of selection - low fitness
    fit_hi = 0.3                        # probability of selection - high fitness
    fitness_states = [fit_lo, fit_hi]
    # this creates the two states

    r0_fitness = fitness_calculation(r0_pop, fitness_states, fitness_function)
    # this calculates the fitness of R0

## Selection

    r1_pop = selection(r0_pop,r0_fitness,trials,sequence_length)
    # this carries out the selection

    noise_level = 0.1
    r0_noise = sel_noise(noise_level, r0_pop, sequence_length)
    # this generates R0 noise
    r1_pop = convert(SharedArray, vcat(r0_noise, r1_pop))
    # this concatenates the matrices

    r1_pop = noise_amp(r1_pop, r1_pop_size)
    # this amplifies the noise

    error_rate = 0.01                    # probability of mutation per position
    mutation!(r1_pop, error_rate, amino_acids)
    # this introduces mutations at random to create noise in the recovered dataset


## Analysing the resulting population (R1)
    r1_fitness = fitness_calculation(r1_pop, fitness_states, fitness_function)
    # this calculates the fitness of R1

    r0_metric = Float64(mean(r0_fitness))
    r1_metric = Float64(mean(r1_fitness))
    improvement = round(log2(r1_metric/r0_metric), digits=2)
    # Log2 chosen here to make the change more obvious

## Showing changes in population
    r1_dist = display(plot([r0_fitness, r1_fitness], bins = 0:0.05:0.75, seriestype = [:barhist :stephist],
    normed=:probability, fillcolor =[:lightgrey :blue], linecolor=[:lightgrey :blue],
    thickness_scaling = 2, label = ["R0" "R1"]))

#= The code in atom is causing race conditions. For instance, improvement and r1_dist
on a "run all" returns 0, while on a step-by-step it returns the correct answer.
It is also incredibly heavy on memory with 1 1e5 scale using up to 34 GiB.
At that scale, in my machine I am getting around 20s running time to this point =#

## Extracting positional frequencies from populations

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

                information_xixjxk_r[i,j,k]    = (entropy_xi + entropy_xj + entropy_xk)
                                - (entropy_xixj + entropy_xixk + entropy_xjxk)
                                + (entropy_xixjxk);
            end
        end
    end

    return information_xixjxk_r
end



## Delta MI in 2D
r0_mi = mi_xy2d(r0_pop)
r1_mi = mi_xy2d(r1_pop)

delta_mi = r1_mi - r0_mi

Heat_map = heatmap(1:size(delta_mi,1),
    1:size(delta_mi,2), delta_mi,
    c=cgrad([:blue, :white, :white, :red]),
    xlabel="aa", ylabel="aa", title="dMI")

## Delta MI in 3D

r0_mi3d = mi_xyz3d(r0_pop)
r1_mi3d = mi_xyz3d(r1_pop)

delta_mi3d = r1_mi3d - r0_mi3d


## Simple analysis of MI in 3D
simple_deltami3d = reshape(delta_mi3d, (length(delta_mi3d)))
simple_deltami3d = sort(simple_deltami3d)

plot([simple_deltami3d], bins = simple_deltami3d[begin]:0.002:simple_deltami3d[end], seriestype = [:barhist], fillcolor =[:orange], linecolor=[:orange],
thickness_scaling = 1, label = "deltaMI 3D") 

findall(x->x==simple_deltami3d[1], delta_mi3d)

## Creating subpopulations and conditional probability analyses

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

    #this function returns all p(x,y) probabilities for a given z residue

    pop_sort = sortslices(population, dims = 1, by= a -> a[z])
    indexes = z_subsets(population, z, amino_acids, sequence_length)

    pxiyi_givenz = zeros(Float64, length(amino_acids), length(amino_acids), length(amino_acids))
    for zi = 1:length(amino_acids)
        pxiyi_givenz[:,:,zi] = freq_xy(convert(SharedArray{Char, 2},pop_sort[indexes[zi]:indexes[zi+1], :]), xy[1],
            xy[2], amino_acids)
    end

    return pxiyi_givenz
end

function delta_pxy_givenz(r1_population::SharedArray{Char, 2}, 
    r0_population::SharedArray{Char, 2}, 
    xy::Array{Int64,1}, z::Int64, 
    amino_acids::Array{Char, 1})

    # this is just a function to minimise errors when looking at the 
    # change in conditional probability P(x,y|z)

    delta_p = zeros(Float64, length(amino_acids),length(amino_acids),length(amino_acids))
    delta_p = pxy_givenz(r1_population, xy ,z, amino_acids) - pxy_givenz(r0_population, xy,z, amino_acids)

    return delta_p
end

xy = [3,4]
z = 5

maps = delta_pxy_givenz(r1_pop, r0_pop, xy, z, amino_acids)
plot(heatmap(maps[:,:,1]), heatmap(maps[:,:,2]), heatmap(maps[:,:,3]), heatmap(maps[:,:,4]),
    heatmap(maps[:,:,5]), heatmap(maps[:,:,6]), heatmap(maps[:,:,7]), clims = (minimum(maps), maximum(maps)), layout = 7)
    # the heatmaps here are limited to the 7 aa currently in the amino_acids - needs to be automated

## Simple analysis of delta p(x,y)|z
simple_deltapxy = reshape(maps, (length(maps)))
simple_deltapxy = sort(simple_deltapxy)
plot([simple_deltapxy], bins = simple_deltapxy[begin]:0.001:simple_deltapxy[end], seriestype = [:barhist], fillcolor =[:lightgreen], linecolor=[:green],
thickness_scaling = 1, label = "delta p(x,y|z)") 

## Simple interrogation of the simple_deltapxy tensor
# it basically recapitulates positive selection functions
combinations = findall(x->x==simple_deltapxy[end], maps)
println("given ", z, " as ", amino_acids[combinations[1][3]], " the following combinations are of interest:")
println("residue ", xy[2], " = ", amino_acids[combinations[1][1]])
println("residue ", xy[1], " = ", amino_acids[combinations[1][2]])


## MI of xy|z - calculated here as the sum of MIs for one Z
#= I don't think this is the fastest way to do this but I am too curious 
to not set it up an see what it looks like. I'll create an array of tensors
and then use one to vizualise the MI and sum of MI =#

function z_subsets(population::SharedArray{Char, 2},
    amino_acids::Array{Char, 1}, 
    sequence_length::Int64)

    # This returns an array containing the indexes for each
    # residue in a sorted array (x for residues, y for aa)

    indexes = zeros(Int64, length(amino_acids)+1, sequence_length)

    for z = 1:sequence_length
        pop_sort = sortslices(population, dims = 1, by= a -> a[z])
        for aa = 1: length(amino_acids)
            indexes[aa, z] = findfirst(b->b == amino_acids[aa], pop_sort[:,z])
        end
    end

    indexes[end, :] = size(pop_sort, 1) * ones(Int64, sequence_length)
    return indexes
end

function z_subsets(population::SharedArray{Char, 2},
    z::Int64,
    amino_acids::Array{Char, 1}, 
    sequence_length::Int64)

    # This returns an array containing the indexes for one
    # residue (z) in a sorted vector

    indexes = zeros(Int64, length(amino_acids)+1)

    pop_sort = sortslices(population, dims = 1, by= a -> a[z])
    for aa = 1: length(amino_acids)
        indexes[aa] = findfirst(b->b == amino_acids[aa], pop_sort[:,z])
    end
    indexes[end] = size(pop_sort, 1)
    return indexes
end

function givenz_mi_xy(population::SharedArray{Char, 2},
    z::Int64,
    sequence_length::Int64,
    amino_acids::Array{Char, 1})

    # this function gives an MI tensor for xy for each amino_acid
    # at a given position z. Thinking about how these values are calculated
    # I reckon a calculation as for the MI functions will be more accurate

    pop_sort = sortslices(population, dims = 1, by= a -> a[z])
    givenz_mixy = zeros(Float64, sequence_length, sequence_length, length(amino_acids))
    pop_index = z_subsets(population, z, amino_acids, sequence_length)

    for aa = 1: length(amino_acids)
        subpop = convert(SharedArray{Char, 2},pop_sort[pop_index[aa]:pop_index[aa+1], :])
        givenz_mixy[:,:,aa] = mi_xy2d(subpop)
    end

    return givenz_mixy
end


# calculates all the required subsets (to save memory)
r0_index = z_subsets(r0_pop, amino_acids, sequence_length)
r1_index = z_subsets(r1_pop, amino_acids, sequence_length)

# calculates dMI(x,y|z)
r0_z1 = givenz_mi_xy(r0_pop,5, sequence_length, amino_acids)
r1_z1 = givenz_mi_xy(r1_pop,5, sequence_length, amino_acids)
deltami_z1 = r0_z1 - r1_z1

maps2 = deltami_z1
plot(heatmap(maps2[:,:,1]), heatmap(maps2[:,:,2]), heatmap(maps2[:,:,3]), heatmap(maps2[:,:,4]),
    heatmap(maps2[:,:,5]), heatmap(maps2[:,:,6]), heatmap(maps2[:,:,7]), clims = (minimum(maps2), maximum(maps2)), layout = 7)
# same constraints to the previous tiled heatmap apply here

# just summing up across the possible amino acids - I don't think
# the maths is correct here.
collapse_mi = zeros(Float64, sequence_length, sequence_length)
for aa1 = 1:sequence_length
    for aa2 = 1:sequence_length
        collapse_mi[aa1, aa2] = sum(deltami_z1[aa1, aa2, :])
    end
end

Heat_map = heatmap(1:size(collapse_mi,1),
    1:size(collapse_mi,2), collapse_mi,
    c=cgrad([:blue, :white, :white, :red]),
    xlabel="aa", ylabel="aa", title="dMI | z for res=5")
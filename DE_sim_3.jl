## Directed evolution simulation and analysis
# July/20 v3.0


# The idea is to create a simulation of directed evolution in which:
# - given a set of parameters, a new population can be created (R0)
# - given a set of stochastic criteria, R0 can be selected
# - The resulting population (R1) includes a component of noise

using BenchmarkTools
using Plots
using SharedArrays
using DelimitedFiles


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
    for position = 1:length(position_composition)
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
    fitness_states::Array{Float64,1},
    fitness_function::Array{Tuple{Array{Int64,1},Array{Array{Char,1},1},Array{Float64,1}},1})
    
    # the function helps create a list of fitness criteria
    fit_comb = (residues, high_fit_combinations, fitness_states)
    return vcat(fitness_function, fit_comb)
end

function fitness_seq(sequence::Array{Char, 1},
    fitness_function::Array{Tuple{Array{Int64,1},Array{Array{Char,1},1},Array{Float64,1}},1})

# this function sums up the fitness of a single sequence

    fitness = 0.0
    for criteria = 1:length(fitness_function)
        if sequence[fitness_function[criteria][1]] in fitness_function[criteria][2]
            fitness += fitness_function[criteria][3][end]
            # Still keeping a 2-state fitness contribution for now
        else
            fitness += fitness_function[criteria][3][1]
        end
    end
    return fitness
end

function fitness_calculation(r0_pop::SharedArray{Char, 2},
    fitness_function::Array{Tuple{Array{Int64,1},Array{Array{Char,1},1},Array{Float64,1}},1})

# this function calculates the individual sequence fitness for a population

    fitness = SharedArray{Float32, 1}(size(r0_pop, 1), init = 0.0)     # Creates an empty fitness matrix
    for entry = 1:size(r0_pop,1)
        sequences = r0_pop[entry, :]
        fitness[entry,1] = fitness_seq(sequences, fitness_function)
    end
    return fitness

end

function fitness_to_recovery(fitness_list::SharedArray{Float32, 1})
    
    # this function converts fitness scores into a probability of recovery
    # equation 1 = logistic equation. For the values used here,
    # scores between -1 and 3 have the biggest impact on selection
    
    selection_p = zeros(Float32, length(fitness_list))
    x0 = 0.5    # point of maximum selection
    L = 0.25  # Maximum recovery 
    k = 3     # how digital the selection is
    
    for a = 1 : length(fitness_list)
        selection_p[a] = L/(1+exp(-k*(fitness_list[a]-x0)))
    end
    return selection_p
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

    for position = 1:length(population)
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
    position_composition = position!(6, ['A','D', 'G', 'F'], position_composition)
    position_composition = pos_criteria_check(position_composition, sequence_length, amino_acids)
    # This sets the positional requirements in R0
    r0_pop = r0_library(position_composition, r0_pop_size, sequence_length)
    # This generates the random R0

    writedlm( "r0_pop.csv",  r0_pop, ',')
    # Commits the starting population to file

## Fitness landscape functions
    fitness_function = Tuple{Array{Int64,1},Array{Array{Char,1},1},Array{Float64,1}}[]
    # Resets the functional landscape

    #fitness_function = fit_criteria!([1], [['C']], [0.0, 1.0], fitness_function)
    fitness_function = fit_criteria!([1], [['A'],['C']], [0, 0.5], fitness_function)
    #fitness_function = fit_criteria!([1, 2, 3], [['C', 'D', 'D']], [0.0, 1.0], fitness_function)
    #fitness_function = fit_criteria!([1, 2, 3], [['C', 'D', 'D'],['A', 'D', 'D'],['A', 'E', 'E'],['A', 'D', 'E']], [-0.1, 0.5], fitness_function)
    #fitness_function = fit_criteria!([3, 4, 5], [['D', 'G', 'G'],['E', 'D', 'D']], [0.0, 2.0], fitness_function)
    #fitness_function = fit_criteria!([3, 4, 5], [['D', 'G', 'G'],['E', 'D', 'D'],['D', 'A', 'A'],['E', 'A', 'A']], [-0.1, 0.5], fitness_function)
    fitness_function = fit_criteria!([3, 4], [['D', 'G'],['E', 'D']], [0.0, 2.0], fitness_function)
    fitness_function = fit_criteria!([2, 3], [['C', 'D'],['C', 'H']], [0.0, -0.3], fitness_function)

## Calculating fitness and probability of selection for all sequences    
    r0_fitness = fitness_calculation(r0_pop, fitness_function)
    # this calculates the fitness of R0
    r0_fitness = convert(SharedArray{Float32,1},fitness_to_recovery(r0_fitness))
    # this converts fitness to probability of selection

    println("Probabilities of recovery are between ", minimum(r0_fitness), " and ", maximum(r0_fitness)) # sanity check 


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

    writedlm( "r1_pop.csv",  r1_pop, ',')
    # Commits the selected population to file

## Analysing the resulting population (R1)
    r1_fitness = fitness_calculation(r1_pop, fitness_function)
    # this calculates the fitness of R1
    r1_fitness = convert(SharedArray{Float32,1},fitness_to_recovery(r1_fitness))
    # this converts fitness to probability of recovery

    r0_metric = Float64(mean(r0_fitness))
    r1_metric = Float64(mean(r1_fitness))
    improvement = round(log2(r1_metric/r0_metric), digits=2)
    # Log2 chosen here to make the change more obvious

## Showing changes in population
    r1_dist = plot([r0_fitness, r1_fitness], bins = 0:0.01:0.3, seriestype = [:barhist :stephist],
    normed=:probability, fillcolor =[:lightgrey :blue], linecolor=[:lightgrey :blue], thickness_scaling = 2,
    label = ["R0" "R1"])
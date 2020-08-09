## Directed evolution simulation analysis
# August/20 v2.1

#= It seems that the most appropriate analysis should be the Information gained
 of the system based on the knowledge of a position.
 This translates (2006 Moore et al) to:
 IG(ABC) = I(A;B|C) - I(A;B)
 I implement that below =#
 
using BenchmarkTools
using Plots
using SharedArrays
using DelimitedFiles
using StatsBase

## Dependencies from selection file:
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H' , 'I', 'K', 'L', 'M']#, 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # Not using a 'packaged' amino acid symbol to allow unnatural amino acids and deletions to be included

## Importing data from selection
    r0_pop = readdlm("r0_pop.csv", ',', Char)
    r0_pop = convert(SharedArray, r0_pop)

    r1_pop = readdlm("r1_pop.csv", ',', Char)
    r1_pop = convert(SharedArray, r1_pop)

    sequence_length = size(r0_pop, 2)

## Functions
function freq_xy_z(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1},
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
    return prob_trio
end

function info_xy_givenz(xyz_distribution::Array{Float64, 3}, amino_acids::Array{Char, 1})

    info = 0
    for a = 1 : length(amino_acids)
        pz   = sum(xyz_distribution[:,:,a])
        pxyz = xyz_distribution[a,a,a]
        pxz  = sum(xyz_distribution[a,:,a]) 
        pyz  = sum(xyz_distribution[:,a,a])
        if pxyz != 0
            info += pxyz * log2((pz*pxyz)/(pxz*pyz))
        end
    end
    return info
end

function info_xy(xyz_distribution::Array{Float64, 3}, amino_acids::Array{Char, 1})
    
    info = 0
    for a = 1 : length(amino_acids)
        pxy   = sum(xyz_distribution[a,a,:])
        px    = sum(xyz_distribution[a,:,:])
        py    = sum(xyz_distribution[:,a,:])
        if pxy != 0
            info += pxy * log2((pxy)/(px*py))
        end
    end
    return info
end

function infogained_xyz(population::SharedArray{Char, 2}, amino_acids::Array{Char, 1})

    ig_3d = zeros(Float64, size(population, 2), size(population, 2),size(population, 2))
    for res1 = 1: size(population, 2)
        for res2 = 1: size(population, 2)
            for res3 = 1: size(population, 2)
                xyz_distribution = freq_xy_z(population, amino_acids, res1, res2, res3);
                xy_givenz = info_xy_givenz(xyz_distribution, amino_acids)
                xy_only = info_xy(xyz_distribution, amino_acids)
                ig_3d[res1,res2,res3] = xy_givenz - xy_only
            end
        end
    end
    return ig_3d
end

function remove_self(info_gained::Array{Float64, 3})
    
    for res1 = 1:size(info_gained, 1)
        for res2 = 1:size(info_gained, 2)
            for res3 = 1:size(info_gained, 3)
                if res1==res2 || res2==res3 || res1==res3
                    info_gained[res1, res2, res3] = 0
                end
            end
        end
    end
    
    return info_gained
end

## Calculations
    infogained = infogained_xyz(r1_pop, amino_acids);
    relevant_info = remove_self(infogained);
    
## Simple analysis of IG

    ig = reshape(relevant_info, (length(relevant_info)))
    ig = sort(ig)
    ## Most positive hits
    println("The following interactions may be of importance:")
    for hit = 1:10
        combinations = findlast(x->x==ig[end+1-hit], relevant_info)
        println("$hit :  IG = I($(combinations[1]);$(combinations[2])|$(combinations[3])) - I($(combinations[1]);$(combinations[2])) -> $(ig[end+1-hit])")
    end
    ## Most negative hits
    println("The following interactions may be of importance:")
    for hit = 1:10
        combinations = findlast(x->x==ig[hit], relevant_info)
        println("$hit :  IG = I($(combinations[1]);$(combinations[2])|$(combinations[3])) - I($(combinations[1]);$(combinations[2])) -> $(ig[hit])")
    end

## Graphic output    
    plot([ig], bins = ig[begin]:0.001:ig[end],
            seriestype = [:barhist], fillcolor =[:orange], linecolor=[:orange],thickness_scaling = 1,
            label = ["IG = I(X;Y|Z) - I(A;B)"]) 

#= That is very similar to some of the implementations carried out in the DE_analysis file. It remains
unclear to me how the relevant data would be extracted from the dataset. Other metrics have given me results
that were more intuitive. =#

using Random, Distributions, SparseArrays, ProgressMeter
include("cycles.jl")

"""
    sample_row_and_col_sums(T::Array{<:Integer})

Given a contingency table `T`, perform the first step of the lumped Burnside process. 
The function samples row and column sums by applying stick breaking to the entries of `T`.  

# Arguments
- `T::Array{<:Integer}`: The input contingency table.

# Returns
- `r::Dict{Integer, Array{Integer}}`: Dictionary of row sums.
- `c::Dict{Integer, Array{Integer}}`: Dictionary of column sums.
- `active::Set{Integer}`: Set of active cycle lengths.
"""
function sample_row_and_col_sums(T::Array{<:Integer})
    I, J = size(T)
    r = Dict{Integer, Array{Integer}}()
    c = Dict{Integer, Array{Integer}}()
    active = Set{Integer}()
    for i in 1:I
        for j in 1:J
            if T[i, j] > 0
                cycles = random_cycles(T[i, j])
                for l in cycles
                    if l in active
                        r_ = r[l]
                        r_[i] += 1
                        r[l] = r_

                        c_ = c[l]
                        c_[j] += 1
                        c[l] = c_
                    else 
                        r_ = zeros(Integer, I)
                        r_[i] = 1
                        r[l] = r_

                        c_ = zeros(Integer, J)
                        c_[j] = 1
                        c[l] = c_

                        push!(active, l)
                    end
                end
            end
        end
    end
    return r, c, active
end

"""
    sample_Fisher_Yates(r::Vector{<:Integer}, c::Vector{<:Integer})

Samples a contingency table with given row sums `r` and column sums `c` 
from the Fisher-Yates distribution. The Fisher-Yates distribution is also
known as the multivariate hypergeometric distribution.

The resulting table will have the same row and column sums as the input vectors.

# Arguments
- `r::Vector{<:Integer}`: Vector of row sums.
- `c::Vector{<:Integer}`: Vector of column sums.

# Returns
- `T::Array{Integer, 2}`: The sampled contingency table.
"""
function sample_Fisher_Yates(r :: Vector{<:Integer}, c :: Vector{<:Integer})
    I = length(r)
    J = length(c)
    n = sum(r)
    c = copy(c)
    r = copy(r)
    @assert n == sum(c) "Row and column sums must be equal"

    T = zeros(Integer, I, J)
    if n == 0
        return T 
    end

    for i in 1:(I-1) 
        for j in 1:(J-1) 
            T[i, j] = rand(Hypergeometric(c[j], sum(c[(j+1):J]), r[i]))
            r[i] -= T[i,j]
            c[j] -= T[i,j]
        end
        T[i, J] = r[i]
        c[J] -= T[i, J]
        r[i] -= T[i, J]
    end 
    T[I, :] = c

    return T 
end

"""
    sample_table(r::Dict{Integer, Array{Integer}}, 
                 c::Dict{Integer, Array{Integer}}, 
                 active::Set{Integer}, 
                 I::Integer, 
                 J::Integer)

This function performs the second step of the lumped Burnside process.
The function produces a contingency table by sampling from the Fisher-Yates 
distribution for each cycle length l. These cycle lengths are then combined to 
produce the new contingency table.

# Arguments
- `r::Dict{Integer, Array{Integer}}`: Dictionary of row sums.
- `c::Dict{Integer, Array{Integer}}`: Dictionary of column sums.
- `active::Set{Integer}`: Set of active cycle lengths.
- `I::Integer`: Number of rows in the resulting table.
- `J::Integer`: Number of columns in the resulting table.

# Returns
- `T::Array{Integer, 2}`: The sampled contingency table.
"""
function sample_table(r ::Dict{Integer, Array{Integer}}, 
                      c ::Dict{Integer, Array{Integer}}, 
                      active :: Set{Integer}, 
                      I :: Integer, 
                      J :: Integer)
    T = zeros(Integer, I, J)
    for l in active
        X = sample_Fisher_Yates(r[l], c[l])
        T .+= l.*X
    end
    return T
end

"""
    lumped_burnside_step(T0::Array{<:Integer})

Performs one step of the lumped Burnside process on the input contingency table `T0`.

# Arguments
- `T0::Array{<:Integer}`: The input contingency table.

# Returns
- `T::Array{Integer, 2}`: The updated contingency table after one step.
"""
function lumped_burnside_step(T0 :: Array{<:Integer})
    I, J = size(T0)
    r, c, active = sample_row_and_col_sums(T0)
    return sample_table(r, c, active, I, J)
end

"""
    lumped_burnside(T0::Array{<:Integer}, reps::Integer)

Performs the lumped Burnside process on the input contingency table `T0` 
for a specified number of repetitions.

# Arguments
- `T0::Array{<:Integer}`: The input contingency table.
- `reps::Integer`: Number of repetitions.

# Returns
- `Tables::Array{Integer, 3}`: Array of contingency tables generated 
  during the process.
"""
function lumped_burnside(T0 :: Array{<:Integer}, reps :: Integer)
    I, J = size(T0)
    Tables = zeros(Integer, reps, I, J)
    Tables[1, :, :] = T0 
    @showprogress for i in 2:reps
        Tables[i, :, :] = lumped_burnside_step(Tables[i-1,:,:])
    end 
    return Tables
end

"""
    lumped_burnside_feature(fun :: Function, T0::Array{<:Integer}, reps::Integer)

Performs the lumped Burnside process on the input contingency table `T0` 
for a specified number of repetitions and computes a feature for each table 
using the provided function `fun`.

# Arguments
- `fun::Function`: Function to compute a feature from a contingency table.
- `T0::Array{<:Integer}`: The input contingency table.
- `reps::Integer`: Number of repetitions.

# Returns
- `features::Array{Float64}`: Array of computed features.
"""
function lumped_burnside_feature(fun :: Function, T0 :: Array{<:Integer},  reps :: Integer)
    features = zeros(reps)
    features[1] = fun(T0)
    @showprogress for i in 2:reps
        T0 = lumped_burnside_step(T0)
        features[i] = fun(T0)
    end
    return features
end

"""
    chi_sq_stat(T::Array{<:Integer})

Computes the chi-squared statistic for the given contingency table `T`.

# Arguments
- `T::Array{<:Integer}`: The input contingency table.

# Returns
- `chi_sq::Float64`: The chi-squared statistic.
"""
function chi_sq_stat(T :: Array{<:Integer})
    r = dropdims(sum(T, dims = 2), dims = 2)
    c = dropdims(sum(T, dims = 1), dims =1)
    E = r * c' ./sum(r)
    return sum((T .- E).^2 ./ E )
end


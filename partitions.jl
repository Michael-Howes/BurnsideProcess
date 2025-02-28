using Random, SparseArrays, Permutations
include("cycles.jl")

"""
    lumped_burnside_step(a :: SparseVector{<:Integer}) 

Perform one step of the lumped Burnside process for partitions. 

# Arguments
- `a :: SparseVector{<:Integer}`: An integer partition in exponential notation.

# Returns
- `b :: SparseVector{<:Integer}`: The sampled integer partition.
"""
function lumped_burnside_step(a :: SparseVector{<:Integer})
    n = length(a)
    b = spzeros(Integer, n)
    active = a.nzind
    for k in active
        cycles = random_cycles(a[k])
        for C in cycles
            if k >= 1
                j = rand(1:k)
                d = gcd(j, k)
                l = trunc(Integer, (k/d)*C)
                b[l] += d 
            else
                b[C] += 1
            end
        end
    end
    return b
end

"""
    transpose(a :: SparseVector{<:Integer})

Computes the transpose of the integer partition a.

# Arguments
- `a :: SparseVector{<:Integer}`: An integer partition in exponential notation.

# Returns
- `b :: SparseVector{<:Integer}`: The transpose of a.
"""
function transpose(a :: SparseVector{<:Integer})
    n = length(a)
    b = spzeros(Integer, n)
    active = a.nzind
    active = sort(active, rev = true)
    k = 0
    L = length(active)
    for l in 1:(L-1)
        m = active[l]
        k += a[m]
        b[k] = m - active[l+1]
    end
    m = active[L]
    k += a[m]
    b[k] = m 
    return(b)
end

"""
    reflected_burnside_step(a :: SparseVector{<:Integer})

Performs one step of the reflected Burnside process for partitions. 
The function first transposes `a` and then performs one step of the
lumped Burnside process.

# Arguments
- `a :: SparseVector{<:Integer}`: An integer partition in exponential notation.

# Returns
- `b :: SparseVector{<:Integer}`: The sampled integer partition.
"""
function reflected_burnside_step(a :: SparseVector{<:Integer})
    b = transpose(a)
    return lumped_burnside_step(b)
end

"""
    lumped_burnside(n :: Integer, reps :: Integer)

Runs the lumped Burnside process on partitions of size `n` for `reps` steps.
The lumped Burnside process is initialized at the partition `1^n` (i.e. the
partition with `n` parts all equal to `1`) 

# Arguments
- `n::Integer`: The size of the integer partitions.
- `reps::Integer`: The number of steps of the Burnside process

# Returns
- `partitions::SparseArray{Integer, 2}`: The sampled integer partitions of shape 
    n by reps.
"""
function lumped_burnside(n :: Integer, reps :: Integer)
    partitions = spzeros(Integer, n, reps)
    a = spzeros(Integer, n)
    a[1] = n
    partitions[1,1] = n
    for i in 2:reps
        a = lumped_burnside_step(a)
        for k in a.nzind
            partitions[k, i] = a[k]
        end
    end 
    return partitions
end 


"""
    reflected_burnside(n :: Integer, reps :: Integer)

Runs the reflected Burnside process on partitions of size `n` for `reps` steps.
The reflected Burnside process is initialized at the partition `1^n` (i.e. the
partition with `n` parts all equal to `1`) 

# Arguments
- `n::Integer`: The size of the integer partitions.
- `reps::Integer`: The number of steps of the Burnside process

# Returns
- `partitions::SparseArray{Integer, 2}`: The sampled integer partitions of shape 
    n by reps.
"""
function reflected_burnside(n :: Integer, reps :: Integer)
    partitions = spzeros(Integer, n, reps)
    a = spzeros(Integer, n)
    a[1] = n
    partitions[1,1] = n
    for i in 2:reps
        a = reflected_burnside_step(a)
        for k in a.nzind
            partitions[k, i] = a[k]
        end
    end 
    return partitions
end 

"""
    reflected_burnside_feature(fun :: Function, n :: Integer, reps :: Integer)

Runs the reflected Burnside process on partitions of size `n` for `reps` steps.
At each step, the function `fun` is applied to the current partition to compute
a feature.

The reflected Burnside process is initialized at the partition `1^n` (i.e. the
partition with `n` parts all equal to `1`) 

# Arguments
- `fun`: Function to compute a feature from an integer partition.
- `n::Integer`: The size of the integer partitions.
- `reps::Integer`: The number of steps of the Burnside process

# Returns
- `features::Array{Float64}`: Array of computed features..
"""
function reflected_burnside_feature(fun :: Function, n :: Integer, reps :: Integer)
    features = zeros(reps)
    a = spzeros(Integer, n)
    a[1] = n
    features[1] = fun(a)
    for i in 2:reps
        a = reflected_burnside_step(a)
        features[i] = fun(a)
    end 
    return features
end 

"""
    lumped_burnside_feature(fun :: Function, n :: Integer, reps :: Integer)

Runs the lumped Burnside process on partitions of size `n` for `reps` steps.
At each step, the function `fun` is applied to the current partition to compute
a feature.

The lumped Burnside process is initialized at the partition `1^n` (i.e. the
partition with `n` parts all equal to `1`) 

# Arguments
- `fun::Function`: Function to compute a feature from an integer partition.
- `n::Integer`: The size of the integer partitions.
- `reps::Integer`: The number of steps of the Burnside process

# Returns
- `features::Array{Float64}`: Array of computed features..
"""
function lumped_burnside_feature(n :: Integer, fun, reps :: Integer)
    features = zeros(reps)
    a = spzeros(Integer, n)
    a[1] = n
    features[1] = fun(a)
    for i in 2:reps
        a = lumped_burnside_step(a)
        features[i] = fun(a)
    end 
    return features
end 

"""
    unlumped_burnside_step(sigma :: Vector{<:Integer})

Performs one step of the commuting graph walk for permutations. 
The function returns a uniformly sampled permutation `tau` that
commuted with `sigma`.

# Arguments
- `sigma :: Vector{<:Integer}`: A permutation in one-line notation.

# Returns
- `tau :: Vector{<:Integer}`: The sampled permutation.
"""
function unlumped_burnside_step(sigma :: Vector{<:Integer})
    tau = zeros(Int, length(sigma))
    Cycles = sort(cycles(Permutation(sigma)), by = length)

    while !isempty(Cycles)
        l = length(Cycles[1])
        a = 0
        while a < length(Cycles) && length(Cycles[a + 1]) == l
            a += 1
        end
        l_cycles = Cycles[1:a]
        shifts = rand(0:(l-1), a)
        eta = randperm(a)
        for i in 1:a
            tau[l_cycles[i]] = l_cycles[eta[i]][((1:l) .+ shifts[i]) .% l .+ 1]
        end
        Cycles = Cycles[a+1:end]
    end

    return tau
end

"""
    permutation_to_partition(sigma :: Vector{<:Integer})

Maps the permutation `sigma` to the integer partition `a` 
corresponding to the cycle type of `sigma`

# Arguments
- `sigma :: Vector{<:Integer}`: A permutation in one-line notation.

# Returns
- `a :: SparseVector{<:Integer}`: The cycle type of `sigma`
"""
function permutation_to_partition(sigma :: Vector{<:Integer})
    n = length(sigma)
    a = spzeros(Integer, n)
    Cycles = cycles(Permutation(sigma))
    for c in Cycles
        a[length(c)] += 1
    end
    return a 
end

"""
    partition_to_permutation(a :: SparseVector{<:Integer})

Maps the partition `a` to a permutation `sigma` with cycle type
`a`. This map is determinisitic. 

# Arguments
- `a :: SparseVector{<:Integer}`: An integer partition in exponential notation

# Returns
- `sigma :: Vector{<:Intger}`: A permutation with cycle type sigma.
"""
function partition_to_permuation(a :: SparseVector{<: Integer})
    n = length(a)
    sigma = zeros(Integer, n)
    i = 1
    for l in a.nzind
        if l == 1
            sigma[i:(i+a[l]-1)] = i:(i+a[l]-1)
            i += a[l]
        else
            for _ in 1:a[l]
                sigma[i:(i+l-2)] = (i+1):(i+l-1)
                sigma[i+l-1] = i
                i += l
            end
        end
    end
    return sigma
end

"""
    reflected_burnside_sample(n :: Integer, steps :: Integer, reps :: Integer)

Samples an integer partition by running the reflected Burnside process for `steps` steps.
This is repeated `reps` times. 

The reflected Burnside process is initialized at the partition `1^n` (i.e. the
partition with `n` parts all equal to `1`).

# Arguments
- `n::Integer`: The size of the integer partitions.
- `steps::Integer`: The number of steps of the reflected Burnside process.
- `reps::Integer`: The number of repitions.

# Returns
- `partitions::SparseArray{Integer, 2}`: The sampled integer partitions of shape 
    n by reps.
"""
function reflected_burnside_sample(n :: Integer, steps :: Integer, reps :: Integer)
    partitions = spzeros(Integer, n, reps)

    for i in 1:reps
        a = spzeros(Integer, n)
        a[1] = n 
        for _ in 1:steps
            a = transpose(a)
            a = lumped_burnside_step(a)
        end 
        for k in a.nzind
            partitions[k,i] = a[k]
        end
    end
    return partitions
end

"""
    num_ls(a :: SparseVector{<:Integer}, l :: Integer)

Returns the number of parts in the integer partition `a`
of size `l`

# Arguments
- `a::SparseVector{<:Intger}`: An integer partition in exponential notation.
- `l::Integer`: The size of the desired parts

# Returns
- `a[l]::Integer`: The number of parts of size `l`
"""
function num_ls(a :: SparseVector{<:Integer}, l :: Integer)
    return a[l]
end 

"""
    num_parts(a :: SparseVector{<:Integer})

Returns the total number of parts in the integer partition `a`.

# Arguments
- `a::SparseVector{<:Intger}`: An integer partition in exponential notation.

# Returns
- `sum(a)`: The number of parts in `a`.
"""
function num_parts(a :: SparseVector{<:Integer})
    return sum(a)
end 

"""
    num_distinct_parts(a :: SparseVector{<:Integer})

Returns the number of parts of distinct length in the integer partition `a`.
For example, the partition 1^n has n parts but only one distinct length.

# Arguments
- `a::SparseVector{<:Intger}`: An integer partition in exponential notation.

# Returns
- `length(a.nzind)`: The number of distinct parts in `a`.
"""
function num_distinct_parts(a :: SparseVector{<:Integer})
    return length(a.nzind)
end 

"""
    largest_part(a :: SparseVector{<:Integer})

Returns the size of the largest part in the integer partition `a`.

# Arguments
- `a::SparseVector{<:Intger}`: An integer partition in exponential notation.

# Returns
- `maximum(a.nzind)`: The size of the largest part in `a`.
"""
function largest_part(a :: SparseVector{<:Integer})
    return maximum(a.nzind)
end 

"""
    tlargest_part(a :: SparseVector{<:Integer}, t :: Integer)

Returns the size of the t^th largest part in the integer partition `a`.
In particular t_largest_part(a, 1) = largest_part(a).

# Arguments
- `a::SparseVector{<:Intger}`: An integer partition in exponential notation.

# Returns
- `L`: The size of the t^th largest part in `a`.
"""
function t_largest_part(a :: SparseVector{<:Integer}, t :: Integer)
    active = a.nzind 
    active = sort(active, rev = true)
    return active[t]
end


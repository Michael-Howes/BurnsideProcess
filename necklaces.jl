using Random, Primes

"""
    sample_fixed_point(j::Int, n::Int, k=2::Int)::Vector{Int}

Given j in the cyclic group of order n, sample a k-ary vector fixed by j.
"""
function sample_fixed_point(j::Int, n::Int, k=2::Int)::Vector{Int}
    @assert 0 ≤ j "j must be greater than 0."
    @assert j < n "j must be strictly less than n."

    h = gcd(j, n)
    ord = n ÷ h

    base = rand(0:(k-1), h)
    x = repeat(base, ord)
    return x
end

"""
    sample_fixed_point(x::Vector{Int})::Int

Given a vector of length n, sample an element of the cyclic group fixing j.
"""
function sample_stabilizer(x::Vector{Int})
    n = length(x)
    for d in sort(divisors(n))
        if group_action(x, d) == x
            if d == n
                return 0
            end
            ord = n ÷ d
            j = rand(0:(ord-1))
            return j * d
        end
    end
    return 0
end

""" 
    group_action(x::Vector{Int}, j::Int)::Vector{Int}

Cyclically shift x by j.
"""
group_action(x::Vector{Int}, j::Int) = [x[(j+1):end]; x[1:j]]

"""
    burnside_proccess(n::Int, reps::Int, k=2::Int)

Run the Burnside process for the cyclic group acting on k-ary strings.

Return
    xs: List of vectors [x_1, x_2,..., x_reps].
    js: Group elements [j_1, j_2,...,j_reps].
"""
function burnside_proccess(n::Int, reps::Int, k=2::Int)
    x = zeros(Int, n)
    xs = [x]
    js = []
    for _ in 1:(reps-1)
        j = sample_stabilizer(x)
        push!(js, j)
        x = sample_fixed_point(j, n, k)
        push!(xs, x)
    end
    j = sample_stabilizer(x)
    push!(js, j)
    return xs, js
end
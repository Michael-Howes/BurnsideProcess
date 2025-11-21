using Random, Primes, LinearAlgebra, SpecialFunctions, LogExpFunctions

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
        if d == n
            return 0
        end
        if is_fixed(x, d, n)
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
    is_fixed(x::Vector{Int}, j::Int, n::Int)::Bool

Check if x is fixed by j.
"""
is_fixed(x::Vector{Int}, j::Int, n::Int) = all(i -> x[i] == x[i+j], 1:(n-j)) && all(i -> x[i] == x[i+j-n], (n-j+1):n)

"""
    burnside_proccess(n::Int, reps::Int, k=2::Int)

Run the Burnside process for the cyclic group acting on k-ary strings.

Return
    xs: List of vectors [x_1, x_2,..., x_reps].
    js: Group elements [j_1, j_2,...,j_reps].
"""
function burnside_process(n::Int, reps::Int, k=2::Int)
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

"""
    transition_kernel(n, k)
"""
function transition_kernel(n, k)
    C = zeros(Float64, (n, n))
    for i in 0:(n-1), j in 0:(n-1)
        gcdin = gcd(i, n)
        gcdijn = gcd(gcdin, j)

        C[i+1, j+1] = 1 / k^gcdin * sum([num_primatives(d, k) * d / n for d in divisors(gcdijn)])
    end
    return C
end

"""
    μ(n)

Mobius function
"""
function μ(n)
    (n != 1) || return 1
    factors = factor(n)
    Set(values(factors)) == Set([1]) || return 0
    return (-1)^length(factors)
end

"""
    π(n, k)

The stationary distribution on G.
"""
function π(n, k)
    p = zeros(n)
    for i in 0:(n-1)
        p[i+1] = gcd(i, n) * log(k)
    end
    return softmax(p)
end

"""
    num_primatives(n, k)

Return the number of primitive sequences of length `n` with `k` colors.
"""
num_primatives(n, k) = sum(μ(n ÷ d) * k^d for d in divisors(n))


"""
    log_num_primatives(n, k)

Return the logarithm of the number of primitive sequences of length `n` with `k` colors.
"""
function log_num_primatives(n, k)
    if n == 1
        return log(k)
    end
    negative_term = logsumexp(d * log(k) for d in divisors(n) if μ(n ÷ d) == -1)
    positive_term = logsumexp(d * log(k) for d in divisors(n) if μ(n ÷ d) == 1)
    logdiff = logsubexp(positive_term, negative_term)
    return logdiff
end

"""
    log_transition_kernel(n, k)
"""
function log_transition_kernel(n, k)
    log_C = zeros(Float64, (n, n))
    for i in 0:(n-1), j in 0:(n-1)
        gcdin = gcd(i, n)
        gcdijn = gcd(gcdin, j)

        log_C[i+1, j+1] = -gcdin * log(k) + logsumexp(log_num_primatives(d, k) + log(d) - log(n) for d in divisors(gcdijn))
    end
    return log_C
end

"""
    log_transition_kernel(n, k)
"""
function log_lumped_transition_kernel(n, k)
    divs = collect(sort(divisors(n)))
    D = length(divs)

    log_C = zeros(Float64, (D, D))
    for j in 1:D
        b = divs[j]
        if b == n
            log_totient = 0
        else
            factors = factor(n ÷ b)
            log_totient = sum(log1p(-1 / prime) + exponent * log(prime) for (prime, exponent) in factors)
        end
        for i in 1:D
            a = divs[i]
            minab = min(a, b)
            log_C[i, j] = log_totient - a * log(k) + logsumexp(log_num_primatives(d, k) + log(d) - log(n) for d in divisors(minab))
        end
    end
    return log_C
end

function lumped_stationary_distribution(n, k)
    divs = sort(divisors(n))
    p = [log(totient(n ÷ d)) + d * log(k) for d in divs]
    return softmax(p)
end

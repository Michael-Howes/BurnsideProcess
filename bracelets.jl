using Random, Primes, LinearAlgebra, SpecialFunctions, LogExpFunctions

"""
    sample_fixed_point(j::Int, r::Int, n::Int, k=2::Int)::Vector{Int}

Given (j, r) in the dihedral group of order 2n, sample a k-ary vector fixed by (j,r).
"""
function sample_fixed_point(j::Int, r::Int, n::Int, k=2::Int)::Vector{Int}
    @assert 0 ≤ j "j must be greater than 0."
    @assert j < n "j must be strictly less than n."
    @assert 0 ≤ r ≤ 1 "r must be 0 or 1."

    if r == 0
        h = gcd(j, n)
        ord = n ÷ h

        base = rand(0:(k-1), h)
        x = repeat(base, ord)
        return x
    else
        x = zeros(Int, n)
        for i in 1:n
            xi = rand(0:(k-1))
            x[i] = xi
            x[mod(j - i, n)+1] = xi
        end
        return x
    end
end

"""
    sample_stabilizer(x::Vector{Int})::Int

Given a vector of length n, sample an element of the cyclic group fixing j.
"""
function sample_stabilizer(x::Vector{Int})
    n = length(x)
    stab = []
    for j in 0:(n-1), r in 0:1
        if is_fixed(x, j, r)
            push!(stab, (j, r))
        end
    end
    return rand(stab)
end

""" 
    group_action(x::Vector{Int}, j::Int)::Vector{Int}

Cyclically shift x by j.
"""
group_action(x::Vector{Int}, j::Int, r::Int) = r == 0 ? [x[(j+1):end]; x[1:j]] : [x[j:-1:1]; x[end:-1:(j+1)]]

""" 
    is_fixed(x::Vector{Int}, j::Int, n::Int)::Bool

Check if x is fixed by j.
"""
is_fixed(x::Vector{Int}, j::Int, r::Int) = x == group_action(x, j, r)

"""
    burnside_proccess(n::Int, reps::Int, k=2::Int)

Run the Burnside process for the cyclic group acting on k-ary strings.

Return
    xs: List of vectors [x_1, x_2,..., x_reps].
    js: Group elements [j_1, j_2,...,j_reps].
"""
function burnside_process(n::Int, reps::Int, k=2::Int, j0=1::Int, r0=0::Int)
    xs = []
    js = [j0]
    rs = [r0]
    j, r = j0, r0
    for _ in 1:(reps-1)
        x = sample_fixed_point(j, r, n, k)
        push!(xs, x)
        j, r = sample_stabilizer(x)
        push!(js, j)
        push!(rs, r)
    end

    x = sample_fixed_point(j, r, n, k)
    push!(xs, x)
    return xs, js, rs
end

# """
#     transition_kernel(n, k)
# """
# function transition_kernel(n, k)
#     C = zeros(Float64, (n, n))
#     for i in 0:(n-1), j in 0:(n-1)
#         gcdin = gcd(i, n)
#         gcdijn = gcd(gcdin, j)

#         C[i+1, j+1] = 1 / k^gcdin * sum([num_primatives(d, k) * d / n for d in divisors(gcdijn)])
#     end
#     return C
# end

# """
#     μ(n)

# Mobius function
# """
# function μ(n)
#     (n != 1) || return 1
#     factors = factor(n)
#     Set(values(factors)) == Set([1]) || return 0
#     return (-1)^length(factors)
# end

"""
    π(n, k)

The stationary distribution on G.
"""
function π(n, k)
    p = zeros(2 * n)
    # Define p on the rotations
    for i in 0:(n-1)
        p[i+1] = gcd(i, n) * log(k)
    end
    # Define p on the reflections
    if n % 2 == 0
        ndiv2 = n ÷ 2
        p[(n+1):2:end] .= ndiv2 * log(k)
        p[(n+2):2:end] .= (ndiv2 + 1) * log(k)
    else
        p[(n+1):end] .= (n + 1) / 2 * log(k)
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

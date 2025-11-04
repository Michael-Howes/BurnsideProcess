using Permutations, Random, SpecialFunctions, LogExpFunctions

"""
    sample_from_stabilizer(sigma::Permutation, p::Integer, k::Integer)

Performs the first step of the Sylow--Burnside process. That is uniformly sample permutations `h` and `g` such that `h` and `g` are in the Sylow p-subgroup of S_{pk}
and `h^{-1}*sigma*g = sigma`.

Note `k` must satisfy `1<= k < p` and `sigma` must be of length `pk`.

# Arguments
- `sigma::Permutation`: Permutation of length `pk`.
- `p::Integer`: Prime number.
- `k::Integer`: Integer satisfying `1 <= k < p`.

# Return
- `h::Permutation` and `g::Permutation`: A unifrom sample from the stabilizer of `sigma`.
- `a::Integer`: The double cosets containing `sigma` has size `p^a`. 
"""
function sample_from_stabilizer(sigma::Permutation, p::Integer, k::Integer)
    @assert length(sigma) == p * k "Permutation length must be pk."
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    g = collect(1:(p*k))
    a = 2 * k
    inv_sigma = inv(sigma)
    for j in 1:k
        start = (j - 1) * p + 1
        ending = j * p
        conj_cycle = sigma.data[start:ending]

        if (is_cycle_in_sylow_subgroup(conj_cycle, p))
            i = rand(0:(p-1))
            g[start:ending] = start .+ (collect(i:(i+p-1)) .% p)
            a -= 1
        end
    end
    g = Permutation(g)
    return sigma * g * inv_sigma, g, a
end

"""
    is_cycle_in_sylow_subgroup(cycle::Array, p::Integer) -> Bool

Return `true` if `cycle` is a cycle in the Sylow subgroup.

The cycles in the subgroup are either fixed points. Or the cycle has length `p` and is (up to a cyclic shift) of the form

    `[start, start + step, start + (2*step % p), ..., start + ((p-1)*step % p)]``

and `start % p == 1`.
"""
function is_cycle_in_sylow_subgroup(cycle::Array, p::Integer)
    if length(cycle) == 1
        return true
    end
    if length(cycle) != p
        return false
    end
    min_location = argmin(cycle)
    cycle = [cycle[min_location:p]; cycle[1:(min_location-1)]]
    if cycle[1] % p != 1
        return false
    end
    start = cycle[1]
    step = cycle[2] - start
    expected = (collect(0:step:(p-1)*step) .% p) .+ start
    return expected == cycle
end

"""
    sample_from_fixed_points(h::Permutation, g::Permutation, p::Integer, k::Integer) -> Permutation

Uniformly sample a permutation `tau` satisfying `inv(h) * tau * g = tau`. 

The method checks that `h` and `g` are both in the Sylow p-subgroup and that `h` and `g` have the same cycle type.

The permutation `tau` is defined separely as a map from the fixed points of `g` to the fixed points of `h` and as a map from the `p` cycles of `g` to the `p`-cycles of `h`.
- On the fixed points, `g` is a uniformly sampled bijection.
- On the `p`-cycles `tau` cyclically shifts each `p`-cyle of `g` and then maps each cycle to a `p`-cycle of `h`.
"""
function sample_from_fixed_points(h::Permutation, g::Permutation, p::Integer, k::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    g_fixed_points = fixed_points(g)
    h_fixed_points = fixed_points(h)
    @assert length(g_fixed_points) == length(h_fixed_points) "h and g must have the same cycle type."

    tau = collect(1:(p*k))
    tau[g_fixed_points] = shuffle(h_fixed_points)

    p_cycles_g = filter((c) -> length(c) == p, cycles(g))
    p_cycles_h = shuffle(filter((c) -> length(c) == p, cycles(h)))
    num_p_cycles = length(p_cycles_g)
    for i in 1:num_p_cycles
        shift = rand(1:p)
        indexes = collect(shift:(shift+p-1)) .% p .+ 1
        input_cycle = p_cycles_g[i][indexes]
        output_cycle = p_cycles_h[i]
        tau[input_cycle] = output_cycle
    end
    return Permutation(tau)
end


"""
    sylow_burnside(p :: Integer, k :: Integer, reps :: Integer)

Runs the Sylow--Burnside process for the symetric group of size `p*k` started at the identity permutation.

# Arguments
- `p::Integer`: Prime number.
- `k::Integer`: Integer satisfying `1 <= p < k`.
- `reps::Integer`: The number of steps of the Sylow--Burnside process.

# Returns
- `permutations::Array{Permutation}`: Samples from the Burnside process. 
    Contain `reps` permutations each of length `p*k`.
- `sizes::Array{Integers}`: The size of the double cosets containing each sampled permutation.
"""
function sylow_burnside(p::Integer, k::Integer, reps::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    sigma = Permutation(p * k)
    sizes = []
    permutations = []
    for _ in 1:(reps-1)
        h, g, a = sample_from_stabilizer(sigma, p, k)
        push!(permutations, sigma)
        push!(sizes, a)
        sigma = sample_from_fixed_points(h, g, p, k)
    end
    push!(permutations, sigma)
    _, _, a = sample_from_stabilizer(sigma, p, k)
    push!(sizes, a)
    return permutations, sizes
end

"""
    log_num_double_cosets(a::Integer, p::Integer, k::Integer)

Computes the logarithm of the number of Sylow p-double cosets of size `p^a` in S_{pk}.

Based on the formula in 

"""
function log_num_double_cosets(a::Integer, p::Integer, k::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    if !(k <= a <= 2 * k)
        return -Inf
    end
    js = collect((2*k-a):k)
    terms = @. logfactorial((k - js) * p) + logfactorial(js) + 2 * (logfactorial(k) - logfactorial(js) - logfactorial(k - js)) + js * (log(p) + log(p - 1)) + logfactorial(js) - logfactorial(2 * k - a) - logfactorial(js - (2 * k - a)) - a * log(p)
    pos_term = logsumexp(terms[1:2:(a-k+1)])
    if a == k
        return pos_term
    else
        neg_term = logsumexp(terms[2:2:(a-k+1)])
        return logsubexp(pos_term, neg_term)
    end
end


"""
    stationary_distribution(p::Integer, k::Integer)

Return the stationary distribution of the Sylow--Burnside process lumped to double coset size.
"""
function stationary_distribution(p::Integer, k::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    log_pi = zeros(2 * k)
    for a in k:2k
        log_pi[a] = log_num_double_cosets(a, p, k)
    end
    return softmax(log_pi)
end

p = 5
k = 1
sylow_burnside(p, k, 10)
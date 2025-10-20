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
    g = Permutation(p * k)
    a = 2 * k
    for j in 1:k
        eta = collect(1:(p*k))
        start = (j - 1) * p + 1
        ending = p * j
        eta[start:ending-1] = collect(start+1:ending)
        eta[ending] = start
        eta = Permutation(eta)
        eta_conj = inv(sigma) * eta * sigma
        is_in_H = is_in_sylow_subgroup(eta_conj, p, k)
        if (is_in_H)
            i = rand(0:(p-1))
            g = g * eta_conj^i
            a -= 1
        end
    end
    return sigma * g * inv(sigma), g, a
end

"""
    is_in_sylow_subgroup(sigma::Permutation, p::Integer, k::Integer) -> Bool

Return `true` if `sigma` (of length `p*k`) lies in the Sylow p-subgroup used by this module.

Preconditions:
- `sigma` must have length `p*k`.
- `k` must be between `1` and `p-1`.

Checks performed:
- every cycle of `sigma` has length either `1` or `p`.
- each `p`-cycle is supported on a block whose first element is congruent to `1 (mod p)`
  and whose entries follow a constant step (mod `p`) so that the cycle is of the
  form {b, b+step, b+2*step, ..., b+(p-1)*step}.

Returns `false` if any condition fails, otherwise `true`.
"""
function is_in_sylow_subgroup(sigma::Permutation, p::Integer, k::Integer)
    @assert length(sigma) == p * k "Permutation length must be pk."
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    for c in cycles(sigma)
        if length(c) != 1 && length(c) != p
            return false
        end
        if length(c) == p
            start = c[1]
            if start % p != 1
                return false
            end
            step = c[2] - c[1]
            expected = (collect(0:step:((p-1)*step)) .% p) .+ start
            if expected != c
                return false
            end
        end
    end
    return true
end

"""
    sample_from_fixed_points(h::Permutation, g::Permutation, p::Integer, k::Integer) -> Permutation

Construct a random permutation `tau` that maps the fixed points of `g` to the fixed points of `h`
and maps each `p`-cycle of `g` bijectively to a `p`-cycle of `h` up to a random cyclic rotation.

Preconditions:
- `1 <= k < p`.
- `h` and `g` must be elements of the Sylow `p`-subgroup.

Behavior:
- Fixed points of `g` are assigned to fixed points of `h` by a uniform shuffle.
- The `p`-cycles of `g` are paired with a uniformly shuffled list of the `p`-cycles of `h`;
  for each pair a random cyclic shift is applied when mapping entries of the `g`-cycle
  to the `h`-cycle.

Returns the resulting permutation `tau`.
"""
function sample_from_fixed_points(h::Permutation, g::Permutation, p::Integer, k::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."

    tau = collect(1:(p*k))
    tau[fixed_points(g)] = shuffle(fixed_points(h))

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
- `k::Integer`: Integer satisfying `1 <= p <k`.
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
        log_num = logexpm1(pos_term - neg_term) + neg_term
        return log_num
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
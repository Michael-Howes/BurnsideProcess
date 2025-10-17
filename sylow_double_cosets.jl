using Permutations, Random

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
        if (is_in_sylow_subgroup(eta_conj, p, k))
            i = rand(0:(p-1))
            g = g * eta^i
            a -= 1
        end
    end
    return sigma * g * inv(sigma), g, a
end

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

function sample_from_fixed_points(h::Permutation, g::Permutation, p::Integer, k::Integer)
    @assert 1 <= k "k must satisfy 1 <= k."
    @assert k < p "k must satisfy k < p."
    @assert is_in_sylow_subgroup(h, p, k) "h must be in the p-Sylow subgroup."
    @assert is_in_sylow_subgroup(g, p, k) "g must be in the p-Sylow subgroup."

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

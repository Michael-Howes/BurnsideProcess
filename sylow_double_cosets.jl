using Permutations

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
"""

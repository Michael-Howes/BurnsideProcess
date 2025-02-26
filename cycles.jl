using Random

"""
    Random_cyles(n::Integer)

Samples the cycle lengths of a uniform permutation of
size `n`. The cycle lengths are sampled using discrete
stick breaking.

# Arguments:
- `n::Integer`: The size of the permutation.

# Returns:
- `cycle_type::Array{::Integer}`: A list of cycle lengths.
"""
function random_cycles(n::Integer)
    @assert n > 0 "n must be positive"
    k = n
    cycle_type = zeros(Integer, 0)
    while k > 0
        cycle = rand(1:k)
        append!(cycle_type, cycle)
        k -= cycle
    end
    return cycle_type
end

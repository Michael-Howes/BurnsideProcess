using Random, Primes

function sample_fixed_point(j::Int, n::Int, k=2::Int)
    @assert 0 ≤ j "j must be greater than 0."
    @assert j < n "j must be strictly less than n."

    h = gcd(j, n)
    ord = n ÷ h

    base = rand(0:(k-1), h)
    x = repeat(base, ord)
    return x
end

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

group_action(x::Vector{Int}, j::Int) = [x[(j+1):end]; x[1:j]]

sample_stabilizer(zeros(Int, 60))
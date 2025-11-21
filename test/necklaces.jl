using Test
include("../necklaces.jl")

@testset "Group action testset" begin
    x = [1, 2, 3, 4]
    j = 2
    @test group_action(x, j) == [3, 4, 1, 2]

    @test group_action(x, 0) == x

    @test group_action(x, 4) == x
end

@testset "Fixed by group action testset" begin
    x = [1, 2, 1, 2]
    @test is_fixed(x, 0, 4)
    @test !is_fixed(x, 1, 4)
    @test is_fixed(x, 2, 4)
    @test !is_fixed(x, 3, 4)


end

@testset "Sample fixed points testset" begin
    n = 10
    j = 5
    x = sample_fixed_point(j, n)
    @test length(x) == n
    @test group_action(x, j) == x

    n = 10
    j = 0
    x = sample_fixed_point(j, n)
    @test length(x) == n
    @test group_action(x, j) == x

    n = 10
    j = 1
    x = sample_fixed_point(j, n)
    @test length(x) == n
    @test group_action(x, j) == x
    @test x == repeat([x[1]], n)
end

@testset "Sample stabilizer testset" begin
    x = repeat([0, 1], 5)
    j = sample_stabilizer(x)
    @test j % 2 == 0
    @test group_action(x, j) == x
    @test j < 10
    @test 0 ≤ j
end

@testset "Burnside process testset" begin
    n = 10
    reps = 50
    xs, js = burnside_process(n, reps)
    for (x, j) in zip(xs, js)
        @test group_action(x, j) == x
    end

    n = 12
    reps = 50
    k = 4
    xs, js = burnside_process(n, reps, k)
    for (x, j) in zip(xs, js)
        @test group_action(x, j) == x
    end
end

@testset "Mobius function testset" begin
    expected = [1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
    actual = [μ(k) for k in 1:10]
    @test actual == expected
end

@testset "Number of primative testset" begin
    n = 5
    k = 3
    expected = k^5 - k
    actual = num_primatives(n, k)
    @test expected == actual

    expected = [2, 2, 6, 12, 30, 54, 126, 240, 504]
    actual = [num_primatives(n, 2) for n in 1:length(expected)]
    @test actual == expected
end

@testset "Transition kernel testset" begin
    n = 8
    k = 3
    P = transition_kernel(n, k)
    @test P * ones(n) ≈ ones(n)

    p = 19
    k = 7
    expected = zeros(p, p)
    expected[2:p, 1:p] .= 1 / p
    expected[1, 2:p] .= 1 / p / k^(p - 1)
    expected[1, 1] = 1 - (p - 1) / p / k^(p - 1)
    actual = transition_kernel(p, k)
    @test actual ≈ expected || !isprime(p)
end

@testset "Stationary distribution testset" begin
    n = 8
    k = 3
    p = π(n, k)
    @test sum(p) ≈ 1

    P = transition_kernel(n, k)
    @test p' * P ≈ p'

    n = 50
    k = 2
    p = π(n, k)
    @test sum(p) ≈ 1

    P = transition_kernel(n, k)
    @test p' * P ≈ p'
end

@testset "log number of primatives testset" begin
    n = 5
    k = 3
    expected = log(k^5 - k)
    actual = log_num_primatives(n, k)
    @test expected ≈ actual

    expected = log.([2, 2, 6, 12, 30, 54, 126, 240, 504])
    actual = [log_num_primatives(n, 2) for n in 1:length(expected)]
    for (a, e) in zip(actual, expected)
        @test a ≈ e
    end
end

@testset "log Transition kernel testset" begin
    n = 8
    k = 3
    P = transition_kernel(n, k)
    log_P = log_transition_kernel(n, k)
    @test P ≈ exp.(log_P)

    p = 19
    k = 7
    expected = zeros(p, p)
    expected[2:p, 1:p] .= 1 / p
    expected[1, 2:p] .= 1 / p / k^(p - 1)
    expected[1, 1] = 1 - (p - 1) / p / k^(p - 1)
    actual = exp.(log_transition_kernel(p, k))
    @test actual ≈ expected || !isprime(p)
end

@testset "log lumped transition kernel testset" begin
    n = 9
    k = 2
    actual = log_lumped_transition_kernel(n, k)

    expected = log.([
        6/9 2/9 1/9
        3/18 10/18 5/18
        1/384 5/576 1139/1152
    ])
    for i in 1:3, j in 1:3
        @test expected[i, j] ≈ actual[i, j]
    end

    n = 2^3
    k = 2
    actual = log_lumped_transition_kernel(n, k)

    expected = log.([
        1/2 1/4 1/8 1/8
        1/4 3/8 3/16 3/16
        1/16 3/32 27/64 27/64
        1/256 3/512 27/1024 987/1024]
    )
    for i in 1:4, j in 1:4
        @test expected[i, j] ≈ actual[i, j]
    end

    n = 2^50
    k = 2
    P = exp.(log_lumped_transition_kernel(n, k))
    @test all(sum(P, dims=2) .≈ 1)
end

@testset "lumped stationary distribution testset" begin
    n = 9
    k = 2
    P = exp.(log_lumped_transition_kernel(n, k))
    p = lumped_stationary_distribution(n, k)
    @test p' * P ≈ p'
    @test sum(p) ≈ 1

    n = 8
    k = 3
    P = exp.(log_lumped_transition_kernel(n, k))
    p = lumped_stationary_distribution(n, k)
    @test p' * P ≈ p'
    @test sum(p) ≈ 1

    n = 2^50
    k = 2
    P = exp.(log_lumped_transition_kernel(n, k))
    p = lumped_stationary_distribution(n, k)
    @test p' * P ≈ p'
    @test sum(p) ≈ 1
end
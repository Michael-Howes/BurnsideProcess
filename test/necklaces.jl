using Test
include("../necklaces.jl")

@testset "Group action testset" begin
    x = [1, 2, 3, 4]
    j = 2
    @test group_action(x, j) == [3, 4, 1, 2]

    @test group_action(x, 0) == x

    @test group_action(x, 4) == x
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
    xs, js = burnside_proccess(n, reps)
    for (x, j) in zip(xs, js)
        @test group_action(x, j) == x
    end

    n = 12
    reps = 50
    k = 4
    xs, js = burnside_proccess(n, reps, k)
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
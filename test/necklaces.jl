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

@testset "Burnside process" begin
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
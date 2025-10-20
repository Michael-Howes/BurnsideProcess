using Test, Random, Permutations

include("../sylow_double_cosets.jl")

@testset "sample_from_stabilizer tests" begin
    p1, k1 = 3, 1
    # Identity permutation.
    sigma1 = Permutation(p1 * k1)
    h1, g1, a1 = sample_from_stabilizer(sigma1, p1, k1)

    @test is_in_sylow_subgroup(h1, p1, k1)
    @test is_in_sylow_subgroup(g1, p1, k1)
    @test inv(h1) * sigma1 * g1 == sigma1
    @test a1 == k1

    p2, k2 = 5, 2
    # Identity permutation.
    sigma2 = Permutation(p2 * k2)
    h2, g2, a2 = sample_from_stabilizer(sigma2, p2, k2)

    @test is_in_sylow_subgroup(h2, p2, k2)
    @test is_in_sylow_subgroup(g2, p2, k2)
    @test inv(h2) * sigma2 * g2 == sigma2
    @test a2 == k2

    p3, k3 = 5, 2
    # Permutation in a large double coset.
    sigma3 = Permutation([1, 10, 2, 9, 3, 8, 4, 7, 5, 6])
    h3, g3, a3 = sample_from_stabilizer(sigma3, p2, k2)
    @test is_in_sylow_subgroup(h3, p3, k3)
    @test is_in_sylow_subgroup(g3, p3, k3)
    @test inv(h3) * sigma3 * g3 == sigma3
    @test a3 == 2 * k3
    @test h3 == Permutation(p3 * k3)
    @test g3 == Permutation(p3 * k3)

    p4, k4 = 5, 2
    # Permutation in a medium sized double coset.
    sigma4 = Permutation([10, 6, 9, 7, 8, 1, 2, 3, 4, 5])
    h4, g4, a4 = sample_from_stabilizer(sigma4, p4, k4)
    @test is_in_sylow_subgroup(h4, p4, k4)
    @test is_in_sylow_subgroup(g4, p4, k4)
    @test inv(h4) * sigma4 * g4 == sigma4
    @test a4 == 2 * k4 - 1
end

@testset "is_in_sylow_subgroup tests" begin
    p, k = 5, 2

    eta1 = Permutation([2, 3, 4, 5, 1, 6, 7, 8, 9, 10])
    eta2 = Permutation([1, 2, 3, 4, 5, 7, 8, 9, 10, 6])


    # check all combinations of powers i,j in 0:4.
    for i in 0:4, j in 0:4
        @test is_in_sylow_subgroup(eta1^i * eta2^j, p, k)
    end

    sigma1 = Permutation([2, 3, 1, 4, 5, 7, 6, 8, 9, 10])
    @test !is_in_sylow_subgroup(sigma1, p, k)

    sigma2 = Permutation([2, 3, 4, 1, 5, 6, 7, 8, 9, 10])
    @test !is_in_sylow_subgroup(sigma2, p, k)
end

@testset "sample_from_fixed_points test" begin
    p, k = 5, 2
    # Test h and g each having one p-cycle.
    h1 = Permutation([2, 3, 4, 5, 1, 6, 7, 8, 9, 10])
    g1 = Permutation([1, 2, 3, 4, 5, 10, 6, 7, 8, 9])
    tau1 = sample_from_fixed_points(h1, g1, p, k)
    @test inv(h1) * tau1 * g1 == tau1

    # Test h and g each having two p-cycles.
    h2 = Permutation([2, 3, 4, 5, 1, 10, 6, 7, 8, 9])
    g2 = Permutation([5, 1, 2, 3, 4, 7, 8, 9, 10, 6])
    tau2 = sample_from_fixed_points(h2, g2, p, k)
    @test inv(h2) * tau2 * g2 == tau2
    _, _, a = sample_from_stabilizer(tau2, p, k)
    @test a == k # tau must be in a small double coset.
end

@testset "sylow_burnside test" begin
    p, k = 7, 6
    reps = p^2

    permutations, sizes = sylow_burnside(p, k, reps)

    @test length(permutations) == reps
    @test length(permutations[1]) == p * k
    @test length(sizes) == reps
    @test sizes[1] == k
    @test all(k .<= sizes .<= 2 * k)
end

@testset "log_num_double_cosets test" begin
    p, k = 11, 1
    n1 = 10 # Number of small double cosets.
    n2 = 329890 # Number of large double cosets.

    @test n1 ≈ exp(log_num_double_cosets(k, p, k))
    @test n2 ≈ exp(log_num_double_cosets(k + 1, p, k))

    p, k = 5, 2
    n2 = 32
    n3 = 64
    n4 = 5792
    @test n2 ≈ exp(log_num_double_cosets(k, p, k))
    @test n3 ≈ exp(log_num_double_cosets(k + 1, p, k))
    @test n4 ≈ exp(log_num_double_cosets(k + 2, p, k))

    p, k, = 11, 4
    log_approx = logfactorial(p * k) - (2 * k) * log(p)
    @test log_approx ≈ log_num_double_cosets(2 * k, p, k)
end

@testset "stationary_distribution test" begin
    p, k = 5, 4
    pi = stationary_distribution(p, k)
    @test length(pi) == 2 * k
    @test all(pi[1:(k-1)] .< 1e-12)
    @test sum(pi) ≈ 1
end
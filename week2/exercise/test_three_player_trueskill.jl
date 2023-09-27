include("trueskill.jl")

using Test

# Test the threeplayer function
@testset "threeplayer function" begin
    skill = Gaussian1DFromMeanAndVariance(25.0, 25.0 * 25.0 / (3.0 * 3.0))
    β = 25.0 / 6.0

    # Compute the posterior skills
    (s1, s2, s3) = threeplayer(skill, skill, skill, β)

    @test mean(s1) ≈ 31.311357968443478 atol = 0.0001
    @test sqrt(variance(s1)) ≈ 6.698818832900481 atol = 0.0001
    @test mean(s2) ≈ 24.999999999999986 atol = 0.0001
    @test sqrt(variance(s2)) ≈ 6.2384699253598885 atol = 0.0001
    @test mean(s3) ≈ 18.688642031556533 atol = 0.0001
    @test sqrt(variance(s3)) ≈ 6.698818832900474 atol = 0.0001
end


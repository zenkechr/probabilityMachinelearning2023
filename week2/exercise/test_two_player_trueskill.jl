include("trueskill.jl")

using Test

# Test the twoplayer function
@testset "twoplayer function" begin
    skill = Gaussian1DFromMeanAndVariance(25.0, 25.0 * 25.0 / (3.0 * 3.0))
    β = 25.0 / 6.0

    # Compute the posterior skills
    (s1, s2) = twoplayer(skill, skill, β)

    @test mean(s1) ≈ 29.205220870033607
    @test sqrt(variance(s1)) ≈ 7.194481348831082 atol = 0.0001
    @test mean(s2) ≈ 20.794779129966404 atol = 0.0001
    @test sqrt(variance(s2)) ≈ 7.194481348831082 atol = 0.0001
end

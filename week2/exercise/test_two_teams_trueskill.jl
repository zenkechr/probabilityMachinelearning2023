include("trueskill.jl")

using Test

# Test the twoteams function
@testset "twoteams function" begin
    skill = Gaussian1DFromMeanAndVariance(25.0, 25.0 * 25.0 / (3.0 * 3.0))
    β = 25.0 / 6.0

    # Compute the posterior skills
    (s11, s12, s21, s22) = twoteams(skill, skill, skill, skill, β)

    @test mean(s11) ≈ mean(s12)
    @test mean(s12) ≈ 27.973540193587954 atol = 0.0001
    @test sqrt(variance(s11)) ≈ sqrt(variance(s11))
    @test sqrt(variance(s11)) ≈ 7.784760957252405 atol = 0.0001
    @test mean(s21) ≈ mean(s22)
    @test mean(s22) ≈ 22.026459806412053 atol = 0.0001
    @test sqrt(variance(s21)) ≈ sqrt(variance(s22))
    @test sqrt(variance(s22)) ≈ 7.784760957252405 atol = 0.0001
end


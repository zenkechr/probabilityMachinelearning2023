# Library for Trueskill algorithm
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

include("factors.jl")

# Computes the TrueSkills for a two-player game
function twoplayer(skill1::Gaussian1D, skill2::Gaussian1D, β)
    bag = DistributionBag(Gaussian1D(0, 0))
    factorList = Vector{Factor}()

    # helper function to add factors to a long list
    function addFactor(f)
        push!(factorList, f)
        return (f)
    end

    # Skill variables for the two players
    s1 = add!(bag)
    s2 = add!(bag)

    # Performance variables for the two players
    p1 = add!(bag)
    p2 = add!(bag)

    # Difference of performances
    d = add!(bag)

    # Gaussian prior for the two player skills
    priorF1 = addFactor(GaussianFactor(skill1, s1, bag))
    priorF2 = addFactor(GaussianFactor(skill2, s2, bag))

    # Gaussian likelihood of performance for the two players
    likelF1 = addFactor(GaussianMeanFactor(β * β, p1, s1, bag))
    likelF2 = addFactor(GaussianMeanFactor(β * β, p2, s2, bag))
    diffF = addFactor(WeightedSumFactor(1, -1, p1, p2, d, bag))
    matchF = addFactor(GreaterThanFactor(0, d, bag))

    # run the message passing schedule
    priorF1.update!(1)
    priorF2.update!(1)
    likelF1.update!(1)
    likelF2.update!(1)
    diffF.update!(3)
    matchF.update!(1)
    diffF.update!(1)
    diffF.update!(2)
    likelF1.update!(2)
    likelF2.update!(2)

    # and now compute the log normalization constant
    println("Z = ", exp(logNormalization(factorList, bag)))

    return (bag[s1], bag[s2])
end

# Computes the TrueSkills for a two-team player game
function twoteams(
    skill11::Gaussian1D,
    skill12::Gaussian1D,
    skill21::Gaussian1D,
    skill22::Gaussian1D,
    β,
)
    bag = DistributionBag(Gaussian1D(0, 0))
    factorList = Vector{Factor}()

    # helper function to add factors to a long list
    function addFactor(f)
        push!(factorList, f)
        return (f)
    end

    # Skill variables for the four players
    s11 = # ADD CORRECT CODE HERE
    s12 = # ADD CORRECT CODE HERE
    s21 = # ADD CORRECT CODE HERE
    s22 = # ADD CORRECT CODE HERE

    # Performance variables for the four players
    p11 = # ADD CORRECT CODE HERE
    p12 = # ADD CORRECT CODE HERE
    p21 = # ADD CORRECT CODE HERE
    p22 = # ADD CORRECT CODE HERE

    # Performance for the two team performances
    t1 = # ADD CORRECT CODE HERE
    t2 = # ADD CORRECT CODE HERE

    # Difference of team performances
    d = # ADD CORRECT CODE HERE

    # Gaussian prior for the four player skills
    priorF11 = # ADD CORRECT CODE HERE
    priorF12 = # ADD CORRECT CODE HERE
    priorF21 = # ADD CORRECT CODE HERE
    priorF22 = # ADD CORRECT CODE HERE

    # Gaussian likelihood of performance for the four players
    likelF11 = # ADD CORRECT CODE HERE
    likelF12 = # ADD CORRECT CODE HERE
    likelF21 = # ADD CORRECT CODE HERE
    likelF22 = # ADD CORRECT CODE HERE

    # Weighted sum factor for the team performances
    teamF1 = # ADD CORRECT CODE HERE
    teamF2 = # ADD CORRECT CODE HERE

    # Difference factor for the team performances
    diffF = # ADD CORRECT CODE HERE

    # Match outcome factor
    matchF = # ADD CORRECT CODE HERE

    # message passing schedule
    # ADD CORRECT CODE HERE

    # and now compute the log normalization constant
    println("Z = ", exp(logNormalization(factorList, bag)))

    return (bag[s11], bag[s12], bag[s21], bag[s22])
end

# Computes the TrueSkills for a three-team player game
function threeplayer(skill1::Gaussian1D, skill2::Gaussian1D, skill3::Gaussian1D, β)
    bag = DistributionBag(Gaussian1D(0, 0))
    factorList = Vector{Factor}()

    # helper function to add factors to a long list
    function addFactor(f)
        push!(factorList, f)
        return (f)
    end

    # Skill variables for the three players
    s1 = # ADD CORRECT CODE HERE
    s2 = # ADD CORRECT CODE HERE
    s3 = # ADD CORRECT CODE HERE

    # Performance variables for the three players
    p1 = # ADD CORRECT CODE HERE
    p2 = # ADD CORRECT CODE HERE
    p3 = # ADD CORRECT CODE HERE

    # Difference of pairwise player performances
    d12 = # ADD CORRECT CODE HERE
    d23 = # ADD CORRECT CODE HERE

    # Gaussian prior for the three player skills
    priorF1 = # ADD CORRECT CODE HERE
    priorF2 = # ADD CORRECT CODE HERE
    priorF3 = # ADD CORRECT CODE HERE

    # Gaussian likelihood of performance for the four players
    likelF1 = # ADD CORRECT CODE HERE
    likelF2 = # ADD CORRECT CODE HERE
    likelF3 = # ADD CORRECT CODE HERE

    # Difference factor for the team performances
    diff12F = # ADD CORRECT CODE HERE
    diff23F = # ADD CORRECT CODE HERE

    # Match outcome factor
    match12F = # ADD CORRECT CODE HERE
    match23F = # ADD CORRECT CODE HERE

    # message passing schedule
    priorF1.update!(1)
    priorF2.update!(1)
    priorF3.update!(1)

    likelF1.update!(1)
    likelF2.update!(1)
    likelF3.update!(1)

    Δ = 1e4
    while (Δ > 1e-6)
        Δ = 0
        Δ = max(Δ, )
        # ADD CORRECT INNER UPDATE SCHEDULE CODE HERE
    end

    diff12F.update!(1)
    diff23F.update!(2)

    likelF1.update!(2)
    likelF2.update!(2)
    likelF3.update!(2)

    # and now compute the log normalization constant
    println("Z = ", exp(logNormalization(factorList, bag)))

    return (bag[s1], bag[s2], bag[s3])
end

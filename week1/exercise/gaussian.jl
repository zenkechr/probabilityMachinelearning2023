# Library for 1D Gaussian messages and distribution
#
# 2023 by Ralf Herbrich
# Hasso-Plattner Institute

"""
Data structure that captures the state of an normalized 1D Gaussian. 
In this represenation, we are storing the precision times mean (tau) and the 
precision (rho). This representations allows for numerically stable products of 
1D-Gaussians.
"""
struct Gaussian1D
    tau::Float64            # the precision mean, tau = μ/σ^2 = μ * rho, is the precision adjusted mean
    rho::Float64            # the precision, rho = 1/σ^2, is the inverse of the variance

    # default constructor checking for precision to be non-negative
    Gaussian1D(tau, rho) =
        (rho < 0) ? error("precision of a Gaussian must be non-negative") :
        new(promote(tau, rho)...)
end
# Initializes a standard Gaussian 
Gaussian1D() = Gaussian1D(0, 1)

"""
    Gaussian1DFromMeanVariance(μ,σ2)

Initializes a Gaussian from mean and variance.
precision mean 'tau' is mean divided by variance
precision 'rho' is 1 over variance
"""
Gaussian1DFromMeanVariance(μ, σ2) = Gaussian1D(μ/σ2, 1/σ2) 

"""
    show(io,g)

Pretty-prints a 1D Gaussian
"""
function Base.show(io::IO, g::Gaussian1D)
    if (g.rho == 0.0)
        print(io, "μ = 0, σ = Inf\n")
    else
        print(io, "μ = ", mean(g), ", σ = ", sqrt(variance(g)), "\n")
    end
end


"""
    mean(g)

Returns the mean of the 1D-Gaussian
```julia-repl
julia> mean(Gaussian1D(1,2))
0.5

julia> mean(Gaussian1DFromMeanVariance(1,2))
1.0
```
"""
mean(g::Gaussian1D) = g.tau / g.rho # ADD CORRECT CODE HERE

"""
variance(g)

Returns the variance of the 1D-Gaussian 
```julia-repl
julia> variance(Gaussian1D(1,2))
0.5

julia> variance(Gaussian1DFromMeanVariance(1,2))
2.0
```
"""
variance(g::Gaussian1D) = 1. / g.rho # ADD CORRECT CODE HERE


"""
    absdiff(g1,g2)

Computes the absolute difference of `g1` and `g2` in terms of tau and rho

difference with mean and variance:
X  = X2 - X1
mean = mean2 - mean1
var = var1 + var2

tau = mean / var --> mean = tau * var = tau / rho 
rho = 1/var --> var = 1/rho

difference with tau and rho:
tau = tau2/rho2 - tau1/rho1
rho = rho1 + rho2

# Examples
```julia-repl
julia> absdiff(Gaussian1D(0,1),Gaussian1D(0,2))
1.0

julia> absdiff(Gaussian1D(0,1),Gaussian1D(0,3))
1.4142135623730951
```
"""
absdiff(g1::Gaussian1D, g2::Gaussian1D) = sqrt(abs(g2.rho-g1.rho)) # ADD CORRECT CODE HERE

# show( Gaussian1D(1,2))

print("absdiff\n")
print( absdiff(Gaussian1D(0,1),Gaussian1D(0,3)), "\n" )

"""
*(g1,g2)

Multiplies two 1D Gaussians together and re-normalizes them
# Examples
```julia-repl
julia> Gaussian1D() * Gaussian1D()
μ = 0.0, σ = 0.7071067811865476
```
"""
function Base.:*(g1::Gaussian1D, g2::Gaussian1D)
  return (Gaussian1D(g1.tau+g2.tau, g1.rho+g2.rho))   # ADD CORRECT CODE HERE
end

print("multiplication\n")
show( Gaussian1D() * Gaussian1D() )

"""
/(g1,g2)

Divides two 1D Gaussians from each other
# Examples
```julia-repl
julia> Gaussian1D(0,1) / Gaussian1D(0,0.5)
μ = 0.0, σ = 1.4142135623730951
```
"""
function Base.:/(g1::Gaussian1D, g2::Gaussian1D)
  return (Gaussian1D(g1.tau-g2.tau, g1.rho-g2.rho))   # ADD CORRECT CODE HERE
end

print("division\n")
show(Gaussian1D(0,1) / Gaussian1D(0,0.5))
"""
logNormProduct(g1,g2)

Computes the log-normalization constant of a multiplication of `g1` and `g2`
# Examples
```julia-repl
julia> logNormProduct(Gaussian1D() * Gaussian1D())
c = 0.28209479177387814
```
"""
function logNormProduct(g1::Gaussian1D, g2::Gaussian1D)
  g1_mean = mean(g1)
  g1_var = variance(g1)
  g2_mean = mean(g2)
  g2_var = variance(g2)
  fac = log(2*pi*(g1_var+g2_var)) + ((g1_mean-g2_mean)^2)/(g1_var+g2_var)
  result = - 0.5 * fac
  # prefac = sqrt(g1.rho*g2.rho/(2*pi*(g1.rho+g2.rho)))
  #expterm = -g1.tau^2/(2*g1.rho) - g2.tau^2/(2*g2.rho) + ((g1.tau+g2.tau)^2)/(2*g1.rho+2*g2.rho)
  #result = prefac * exp(expterm)
  return result
end

"""
    logNormRatio(g1,g2)

Computes the log-normalization constant of a division of `g1` with `g2`
# Examples
```julia-repl
julia> logNormRatio(Gaussian1D(0,1) / Gaussian1D(0,0.5))
5.013256549262001
```
"""
function logNormRatio(g1::Gaussian1D, g2::Gaussian1D)
    # ADD CORRECT CODE HERE
    g1_mean = mean(g1)
    g1_var = variance(g1)
    g2_mean = mean(g2)
    g2_var = variance(g2)
    epsilon = eps()
    
    logterm = -0.5*log((g2_var-g1_var+epsilon)/(2*pi*g2_var^2))
    secterm = 0.5 * ((g1_mean-g2_mean)^2)/(g2_var-g1_var+epsilon)
    if isapprox(g2_var-g1_var, 0.0; atol=eps(Float64), rtol=0)
      result = 0.0
    else
      result = logterm + secterm
    end
    
    return result
end


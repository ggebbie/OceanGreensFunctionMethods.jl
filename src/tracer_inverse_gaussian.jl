"""
    TracerInverseGaussian(Γ,Δ)
using LinearAlgebra: NumberArray

The *tracer inverse Gaussian distribution* with mean `Γ` and width `Δ` has probability density function

```math
G(𝐱, \\tau) = \\sqrt{\\frac{\\Gamma^3 }{4 \\pi \\Delta^2 \\tau^3 }} \\exp \\left( - \\frac{\\Gamma (\\tau - \\Gamma)^2}{4 \\Delta ^2 \\tau}\\right) 
```

```julia
TracerInverseGaussian()              # Tracer Inverse Gaussian distribution with unit mean and unit width, i.e. TracerInverseGaussian(1, 1)
TracerInverseGaussian(Γ, Δ)          # Tracer Inverse Gaussian distribution with mean Γ and width Δ

params(d)           # Get the parameters, i.e. (Γ, Δ)
mean(d)             # Get the mean parameter, i.e. Γ
shape(d)            # Get the shape parameter, i.e. Δ
```

External links

* [Compare to Inverse Gaussian distribution on Wikipedia](http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution)
"""
struct TracerInverseGaussian{T<:Number} <: ContinuousUnivariateDistribution
    Γ::T
    Δ::T
    TracerInverseGaussian{T}(Γ::T, Δ::T) where {T<:Number} = new{T}(Γ, Δ)
end

function TracerInverseGaussian(Γ::T, Δ::T; check_args::Bool=true) where {T<:Number}
    @check_args TracerInverseGaussian (Γ, Γ > zero(Γ)) (Δ, Δ > zero(Δ))
    return TracerInverseGaussian{T}(Γ, Δ)
end

TracerInverseGaussian(Γ::Number, Δ::Number; check_args::Bool=true) = TracerInverseGaussian(promote(Γ, Δ)...; check_args=check_args)
TracerInverseGaussian(Γ::Integer, Δ::Integer; check_args::Bool=true) = TracerInverseGaussian(float(Γ), float(Δ); check_args=check_args)
TracerInverseGaussian(Γ::Number; check_args::Bool=true) = TracerInverseGaussian(Γ, one(Γ); check_args=check_args)
TracerInverseGaussian() = TracerInverseGaussian{Float64}(1.0, 1.0)

@distr_support TracerInverseGaussian 0.0 Inf

#### Conversions

function convert(::Type{TracerInverseGaussian{T}}, Γ::S, Δ::S) where {T <: Real, S <: Real}
    TracerInverseGaussian(T(Γ), T(Δ))
end
function Base.convert(::Type{TracerInverseGaussian{T}}, d::TracerInverseGaussian) where {T<:Real}
    TracerInverseGaussian{T}(T(d.Γ), T(d.Δ))
end
Base.convert(::Type{TracerInverseGaussian{T}}, d::TracerInverseGaussian{T}) where {T<:Real} = d

#### Parameters

shape(d::TracerInverseGaussian) = d.Γ^3/(2*d.Δ^2) 
width(d::TracerInverseGaussian) = d.Δ
params(d::TracerInverseGaussian) = (d.Γ, d.Δ)
partype(::TracerInverseGaussian{T}) where {T} = T

# constructor for original Inverse Gaussian
InverseGaussian(d::TracerInverseGaussian) = InverseGaussian(ustrip.(d.Γ), ustrip(shape(d)))

# #### Statistics

mean(d::TracerInverseGaussian) = d.Γ

var(d::TracerInverseGaussian) = d.Γ^3 / shape(d)

skewness(d::TracerInverseGaussian) = 3sqrt(d.Γ / shape(d))

# kurtosis(d::TracerInverseGaussian) = 15d.Γ / d.Δ

# function mode(d::TracerInverseGaussian)
#     Γ, Δ = params(d)
#     r = Γ / Δ
#     Γ * (sqrt(1 + (3r/2)^2) - (3r/2))
# end


# #### Evaluation

function pdf(d::TracerInverseGaussian{T}, x::Number) where T<:Number
    unt = unit(d.Γ) 
    dig = InverseGaussian(d)
    return pdf(dig,ustrip(x))./unt
end

# function logpdf(d::TracerInverseGaussian{T}, x::Real) where T<:Real
#     if x > 0
#         Γ, Δ = params(d)
#         return (log(Δ) - (log2π + 3log(x)) - Δ * (x - Γ)^2 / (Γ^2 * x))/2
#     else
#         return -T(Inf)
#     end
# end

# function cdf(d::TracerInverseGaussian, x::Real)
#     Γ, Δ = params(d)
#     y = max(x, 0)
#     u = sqrt(Δ / y)
#     v = y / Γ
#     z = normcdf(u * (v - 1)) + exp(2Δ / Γ) * normcdf(-u * (v + 1))

#     # otherwise `NaN` is returned for `+Inf`
#     return isinf(x) && x > 0 ? one(z) : z
# end

# function ccdf(d::TracerInverseGaussian, x::Real)
#     Γ, λ = params(d)
#     y = max(x, 0)
#     u = sqrt(λ / y)
#     v = y / Γ
#     z = normccdf(u * (v - 1)) - exp(2λ / Γ) * normcdf(-u * (v + 1))

#     # otherwise `NaN` is returned for `+Inf`
#     return isinf(x) && x > 0 ? zero(z) : z
# end

# function logcdf(d::TracerInverseGaussian, x::Real)
#     Γ, λ = params(d)
#     y = max(x, 0)
#     u = sqrt(λ / y)
#     v = y / Γ

#     a = normlogcdf(u * (v - 1))
#     b = 2λ / Γ + normlogcdf(-u * (v + 1))
#     z = logaddexp(a, b)

#     # otherwise `NaN` is returned for `+Inf`
#     return isinf(x) && x > 0 ? zero(z) : z
# end

# function logccdf(d::TracerInverseGaussian, x::Real)
#     Γ, λ = params(d)
#     y = max(x, 0)
#     u = sqrt(λ / y)
#     v = y / Γ

#     a = normlogccdf(u * (v - 1))
#     b = 2λ / Γ + normlogcdf(-u * (v + 1))
#     z = logsubexp(a, b)

#     # otherwise `NaN` is returned for `+Inf`
#     return isinf(x) && x > 0 ? oftype(z, -Inf) : z
# end

# @quantile_newton TracerInverseGaussian

# #### Sampling

# # rand method from:
# #   John R. Michael, William R. Schucany and Roy W. Haas (1976)
# #   Generating Random Variates Using Transformations with Multiple Roots
# #   The American Statistician , Vol. 30, No. 2, pp. 88-90
# function rand(rng::AbstractRNG, d::TracerInverseGaussian)
#     Γ, λ = params(d)
#     z = randn(rng)
#     v = z * z
#     w = Γ * v
#     x1 = Γ + Γ / (2λ) * (w - sqrt(w * (4λ + w)))
#     p1 = Γ / (Γ + x1)
#     u = rand(rng)
#     u >= p1 ? Γ^2 / x1 : x1
# end

# #### Fit model

# """
# Sufficient statistics for `TracerInverseGaussian`, containing the weighted
# sum of observations, the weighted sum of inverse points and sum of weights.
# """
# struct TracerInverseGaussianStats <: SufficientStats
#     sx::Float64      # (weighted) sum of x
#     sinvx::Float64   # (weighted) sum of 1/x
#     sw::Float64      # sum of sample weight
# end

# function suffstats(::Type{<:TracerInverseGaussian}, x::AbstractVector{<:Real})
#     sx = sum(x)
#     sinvx = sum(inv, x)
#     TracerInverseGaussianStats(sx, sinvx, length(x))
# end

# function suffstats(::Type{<:TracerInverseGaussian}, x::AbstractVector{<:Real}, w::AbstractVector{<:Real})
#     n = length(x)
#     if length(w) != n
#         throw(DimensionMismatch("Inconsistent argument dimensions."))
#     end
#     T = promote_type(eltype(x), eltype(w))
#     sx = zero(T)
#     sinvx = zero(T)
#     sw = zero(T)
#     @inbounds @simd for i in eachindex(x)
#         sx += w[i]*x[i]
#         sinvx += w[i]/x[i]
#         sw += w[i]
#     end
#     TracerInverseGaussianStats(sx, sinvx, sw)
# end

# function fit_mle(::Type{<:TracerInverseGaussian}, ss::TracerInverseGaussianStats)
#     mu = ss.sx / ss.sw
#     invlambda = ss.sinvx / ss.sw  -  inv(mu)
#     TracerInverseGaussian(mu, inv(invlambda))
# end

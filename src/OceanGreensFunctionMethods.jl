module OceanGreensFunctionMethods

using Distributions
using Distributions: @check_args
using Distributions: @distr_support
using DynamicQuantities

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params

#Distributions.InverseGaussian(Γ::Quantity, Δ::Quantity) = InverseGaussian( 
export TracerInverseGaussian
export width
export # re-export from Distributions
    mean, median, quantile, std, var, cov, cor, shape, params

include("tracer_inverse_gaussian.jl")
    
end

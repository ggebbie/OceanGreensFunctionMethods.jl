module OceanGreensFunctionMethods

using Distributions
using Distributions: @check_args
using Distributions: @distr_support
using DimensionalData
using Unitful

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params
import Base: +

export Fluxes
export TracerInverseGaussian
export width
export abyssal_overturning
export intermediate_overturning
export vertical_diffusion
export tracer_flux, tracer_flux_convergence
export projectdir, datadir, srcdir
export # re-export from Distributions
    mean, median, quantile, std, var, cov, cor, shape, params
export # re-export from Base
    +
    
projectdir() = dirname(Base.active_project())
projectdir(args...) = joinpath(projectdir(), args...)

datadir() = joinpath(projectdir(),"data")
datadir(args...) = joinpath(datadir(), args...)

srcdir() = joinpath(projectdir(),"src")
srcdir(args...) = joinpath(srcdir(), args...)

#include("config_units.jl")

include("tracer_inverse_gaussian.jl")

include("pedagogical_tracer_box_model.jl")

end

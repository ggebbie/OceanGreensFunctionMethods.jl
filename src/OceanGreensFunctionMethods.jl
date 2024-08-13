module OceanGreensFunctionMethods

using Distributions
using Distributions: @check_args
using Distributions: @distr_support
using DimensionalData
using DimensionalData: @dim 
using Unitful
using MultipliableDimArrays
using LinearAlgebra
using Downloads
using MAT

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params
import Base: +
import DimensionalData: dims
import LinearAlgebra: eigen

export Fluxes
export TracerInverseGaussian
export width
export abyssal_overturning
export intermediate_overturning
export vertical_diffusion
export advective_diffusive_flux
export convergence, mass_convergence
export tracer_tendency
export projectdir, datadir, srcdir
export meridional_names, vertical_names
export boundary_flux, local_boundary_flux
export linear_probe, mass, uniform
export maximum_timescale, mean_age
export ttd_width, normalized_exponential_decay
export read_tracer_histories
export greens_function
export # re-export from Distributions
    mean, median, quantile, std, var, cov, cor, shape, params
export watermass_fraction
export # re-export from Base
    +
export # re-export from DimensionalData
    dims
export # re-export from LinearAlgebra
    eigen

@dim Eigenmode "eigenmode"

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

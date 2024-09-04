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
using Interpolations
using QuadGK

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params, pdf, InverseGaussian
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
export model_dimensions, boundary_dimensions
export linear_probe, mass #, uniform
export maximum_timescale, mean_age
export adjoint_mean_age
export ttd_width, normalized_exponential_decay
export adjoint_ttd_width
export residence_time
export read_tracer_histories, tracer_source_history
export greens_function
export forward_boundary_propagator
export adjoint_boundary_propagator
export global_ttd, adjoint_global_ttd
export watermass_fraction
export tracer_units
export integrate_forcing
export evolve_concentration
export timestep_initial_condition
export Tracer, Meridional, Vertical, Global 
export transient_tracer_timeseries
export # re-export from Distributions
    mean, median, quantile, std, var, cov, cor, shape, params
export # re-export from Distributions
    InverseGaussian, pdf
export # re-export from Base
    +
export # re-export from DimensionalData
    dims
export # re-export from LinearAlgebra
    eigen

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

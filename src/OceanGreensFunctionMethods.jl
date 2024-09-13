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
using CSV

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params, pdf, InverseGaussian
import Base: +, alignment
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
export ttd_width, normalized_exponential_decay
export residence_time
export read_transient_tracer_histories, read_iodine129_history
export tracer_source_history
export greens_function
export boundary_propagator
export global_ttd
export watermass_fraction
export tracer_units
export integrate_forcing
export evolve_concentration
export timestep_initial_condition
export path_density
export Tracer, Meridional, Vertical, Global 
export tracer_timeseries
export ideal_age
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

function projectdir()
    localproject = dirname(Base.active_project())
    if localproject[end-8:end] == "notebooks"
       return localproject[begin:end-10]
    else
       return(localproject)
    end
end

projectdir(args...) = joinpath(projectdir(), args...)

datadir() = joinpath(projectdir(),"data")
datadir(args...) = joinpath(datadir(), args...)

notebookdatadir() = joinpath(projectdir(),"notebook","data")
notebookdatadir(args...) = joinpath(notebookdatadir(), args...)

srcdir() = joinpath(projectdir(),"src")
srcdir(args...) = joinpath(srcdir(), args...)

#include("config_units.jl")

# credit to Rafael Schouten for properly aligning Unitful quantities 
function Base.alignment(io::IO, x::Unitful.Quantity)
    total = Base.alignment_from_show(io, x)[1]
    parent = Base.alignment(io, ustrip(x))
    return parent[1], parent[2] + total - sum(parent)
end

include("tracer_inverse_gaussian.jl")

include("pedagogical_tracer_box_model.jl")

include("greens_functions.jl")

end

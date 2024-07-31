module OceanGreensFunctionMethods

using Distributions
using Distributions: @check_args
using Distributions: @distr_support
using DimensionalData
using Unitful

import Distributions: mean, median, quantile, std, var, cov, cor, shape, params

export TracerInverseGaussian
export width
export abyssal_overturning
export intermediate_overturning
export vertical_diffusion
export # re-export from Distributions
    mean, median, quantile, std, var, cov, cor, shape, params

include("config_units.jl")

include("tracer_inverse_gaussian.jl")

include("pedagogical_tracer_box_model.jl")

end

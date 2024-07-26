# OceanGreensFunctionMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/dev/)
[![Build Status](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl)

Julia package to complement "A Review of Green's Function Methods in Ocean Circulation Models," by Haine et al. This package goes toward one of the stated goals of the manuscript, namely to make Green's Function methods accessible for learning purposes. Here, we also aim to make a Julia package that is useful and computationally efficient for research purposes.

This package is in an early state and breaking changes are expected.

# Usage

See the tests in `test/runtests.jl` for a more-detailed example usage and a comparison to the existing `InverseGaussian` distribution in the `Distributions.jl` standard package.
```julia
    Γ = 20.0 # mean
    Δ = 20.0 # width

    G = TracerInverseGaussian(Γ, Δ)
	shape(G)
	params(G)
	partype(G)
	mean(G)
	var(G)
```

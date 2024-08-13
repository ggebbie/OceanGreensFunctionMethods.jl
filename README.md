# OceanGreensFunctionMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/dev/)
[![Build Status](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl)

Julia package to complement "A Review of Green's Function Methods in Ocean Circulation Models," by Haine et al. This package has two goal.

1. One of the stated goals of the manuscript is to make Green's Function methods accessible for learning purposes. A computational notebook is provided with this package and detailed below.

2. Here, we also aim to make a Julia package that is useful and computationally efficient for research purposes. The codes are contained in the source directory (`src`) and can be used and imported by other Julia projects. 

This package is in an early state and breaking changes are expected.

# Usage

## Computational Notebook

A Pluto notebook is included (`notebooks/pedagogical_box_model.jl`). Suggested steps to run the notebook (skip steps that you have previously completed):

1. Install julia. I recommend [`juliaup`](https://github.com/JuliaLang/juliaup). To install `juliaup` on a Linux machine, it's as easy as `curl -fsSL https://install.julialang.org | sh`. 

`juliaup status` shows you which Julia versions you have installed and which one is configured as the default.

2. Navigate to the notebook and open julia in the terminal or your favorite editor/IDE. Set up the computing environment.

```julia
Pkg.activate(".") # activate the notebooks reproducible environment
Pkg.instantiate() # download necessary packages if this has not been done before
Pkg.status() # check the status of environment
```
These steps can be accomplished in the REPL package manager mode as well.

3. Start [Pluto](https://plutojl.org/).

```julia
using Pluto
Pluto.run()
```

4. Select the `pedagogical_box_model.jl` file. Exit safe mode and allow the notebook to be run.


## Julia package

This package can be included in other projects by referencing the GitHub URL.
```julia
Pkg.add(url="https://github.com/ggebbie/OceanGreensFunctionMethods.jl")
```

Try out the test suite by switching to the package manager mode at the REPL by typing `]`, then input `test OceanGreensFunctionMethods`.

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

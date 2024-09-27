# OceanGreensFunctionMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/OceanGreensFunctionMethods.jl/dev/)
[![Build Status](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/OceanGreensFunctionMethods.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/OceanGreensFunctionMethods.jl)

Julia package to complement "A Review of Green's Function Methods in Ocean Circulation Models," by Haine et al. This package has two goals.

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
Pkg.activate(".") # activate the notebook's reproducible environment
Pkg.instantiate() # download necessary packages if this has not been done before
Pkg.status() # check the status of environment
```
These steps can be accomplished in the REPL package manager mode as well.

3. Start [Pluto](https://plutojl.org/).

```julia
using Pluto
Pluto.run()
```
 
The notebook should open in a browser or manually enter the URL from the command line into a browser.

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

# Features

1. Self-documented data formats.

Box models and other gridded datasets can be stored in concise arrays that are organized according to physical space or some other organizational system. Geo-scientists are familiar with storing output data to file with meta-data in a self-describing format. Here we use `DimensionalData.jl` to keep data organized during the computational phase of the work as well.

2. Minimizing bookkeeping of obscure indices.

N-dimensional gridded data is naturally stored in N-dimensional arrays, but N-dimensional arrays are not in the right format to perform linear algebra operations with matrices and vectors. The Julia package `MultipliableDimArrays.jl` does the right thing to permit these linear algebra operations and then returns the output in the same human-readable output that the investigator originally provided. Cryptic references to boxes by a given sequential number are minimized in the code, and instead values can be looked up from more easily-interpretable names.

3. Embedding physical units with numerical quantities.

Green's function methods with physical data are a type of Dimensional Analysis as laid out by George W. Hart ("Multidimensional Analysis: Algebras and Systems for Science and Engineering" (Springer-Verlag, 1995). His approach fits nicely into Julia's type system and multiple dispatch using the `Unitful.jl` package. Units do not cause much computational overhead, they remind us of the physical problem at hand, and they provide a means of checking when our equations are not implemented properly.



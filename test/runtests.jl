using Revise
using OceanGreensFunctionMethods
using Distributions
using DynamicQuantities
# using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomadinSets, Unitful
# using DrWatson
# using Plots
using Test

#include("config_units.jl")

@testset "OceanGreensFunctionMethods.jl" begin

    # compare a Tracer Inverse Gaussian distribution
    # with geophysical arguments
    # to an Inverse Gaussian with statistical arguments
    Γ = 20.0 # mean
    Δ = 20.0 # width

    G = TracerInverseGaussian(Γ, Δ)

    # should be identical to this Inverse Gaussian
    G2 = InverseGaussian(mean(G),shape(G))

    @test isequal(shape(G),shape(G2))
    @test !isequal(params(G),params(G2))
    @test isequal(partype(G), partype(G2))
    @test isequal(mean(G), mean(G2))
    @test isequal(var(G), var(G2))

end

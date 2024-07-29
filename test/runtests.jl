using Revise
using OceanGreensFunctionMethods
using Distributions
using DimensionalData
using DimensionalData: @dim
using Unitful
# using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomadinSets, Unitful
# using DrWatson
# using Plots
using Test

#include("config_units.jl")

@testset "OceanGreensFunctionMethods.jl" begin

    @testset "tracer inverse Gaussian distribution" begin
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

    @testset "pedagogical tracer box model" begin

        ENV["UNITFUL_FANCY_EXPONENTS"] = true
        yr = u"yr" # more-readable shorthand for u macro
        Sv = u"Sv" # shorthand for units of sverdrups
        m = u"m"
        km = u"km"

        @dim MeridionalLocation "meridional location"
        @dim VerticalLocation "vertical location"

        # define grid
        meridional_locs = ["High latitudes", "Mid-latitudes", "Low latitudes"]
        vertical_locs = [:Thermocline, :Deep, :Abyssal]
        Ny = length(meridional_locs); Nz = length(vertical_locs)
        
        V_uniform = 1e16m^3 |> km^3 # uniform value of volume for all boxes
        V = DimArray(fill(V_uniform, Ny, Nz),
            (MeridionalLocation(meridional_locs),VerticalLocation(vertical_locs)))

    end
    
end

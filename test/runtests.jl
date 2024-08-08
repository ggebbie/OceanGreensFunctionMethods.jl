using Revise
using OceanGreensFunctionMethods
using Distributions
using DimensionalData
using DimensionalData: @dim
using Unitful
using MultipliableDimArrays
using LinearAlgebra
# using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomadinSets, Unitful
# using DrWatson
# using Plots
using Test

include(srcdir("config_units.jl"))

@testset "OceanGreensFunctionMethods.jl" begin

    @testset "tracer inverse Gaussian distribution" begin
        # compare a Tracer Inverse Gaussian distribution
        # with geophysical arguments
        # to an Inverse Gaussian with statistical arguments
        Γ = 20.0 # mean
        Δ = 20.0 # width

        G = TracerInverseGaussian(Γ, Δ)

        # the Inverse Gaussian distribution in typical statistical toolboxes
        # has different input arguments
        μ = Γ # the same, but different notation
        λ = shape(G) # the shape parameter is the second argument 
        
        G2 = InverseGaussian(μ, λ)

        # both Inverse Gaussians should be the same
        @test isequal(shape(G),shape(G2))
        @test !isequal(params(G),params(G2))
        @test isequal(partype(G), partype(G2)) # parameter type
        @test isequal(mean(G), mean(G2))
        @test isequal(var(G), var(G2))
    end

    @testset "pedagogical tracer box model" begin

        @dim Meridional "meridional location"
        @dim Vertical "vertical location"

        # define grid
        meridional_locs = ["1 High latitudes", "2 Mid-latitudes", "3 Low latitudes"]
        vertical_locs = ["1 Thermocline", "2 Deep", "3 Abyssal"]
        Ny = length(meridional_locs); Nz = length(vertical_locs)
        
        Vol_uniform = 1e16m^3 |> km^3 # uniform value of volume for all boxes

        model_dims = (Meridional(meridional_locs),Vertical(vertical_locs))
        Vol = DimArray(fill(Vol_uniform, Ny, Nz), model_dims)

        Ψ_abyssal = 20Sv
        Fv_abyssal = abyssal_overturning(Ψ_abyssal, model_dims) # volume fluxes

        Ψ_intermediate = 10Sv
        Fv_intermediate = intermediate_overturning(Ψ_intermediate, model_dims) # volume fluxes

        Fv_exchange = 5Sv
        Fv_diffusion = vertical_diffusion(Fv_exchange, model_dims) # volume fluxes

        Fv = Fv_abyssal + Fv_intermediate + Fv_diffusion

        C = ones(model_dims)

        J = advective_diffusive_flux(C, Fv)
        deldotJ = convergence(J)
        deldotFm = mass_convergence(Fv)

        @test all(isapprox.(
            mass_convergence(Fv),
            0.0Tg/s,
            atol=1e-8Tg/s)
        ) # check mass conservation
        
        @test all(isapprox.(
            mass_convergence(Fv_diffusion),
            0.0Tg/s,
            atol=1e-8Tg/s)
        ) # check mass conservation

        @test all(isapprox.(
            mass_convergence(Fv_abyssal),
            0.0Tg/s,
            atol=1e-8Tg/s)
        ) # check mass conservation

        @test all(isapprox.(
            mass_convergence(Fv_intermediate),
            0.0Tg/s,
            atol=1e-8Tg/s)
        ) # check mass conservation

        # other interesting functions
        J_abyssal = advective_diffusive_flux(C, Fv_abyssal)
        deldotJ_abyssal = convergence(J_abyssal)
        
        J_intermediate = advective_diffusive_flux(C, Fv_intermediate)
        deldotJ_intermediate = convergence(J_intermediate)

        @testset "surface boundary" begin 

            # boundary exchange: define the locations affected by boundary fluxes
            meridional_boundary = ["1 High latitudes", "2 Mid-latitudes"]
            vertical_boundary = ["1 Thermocline"]
            boundary_dims = (Meridional(meridional_boundary), Vertical(vertical_boundary))

            Fb = DimArray(hcat([10Sv, 5Sv]), boundary_dims) # boundary flux
            f = ones(boundary_dims) # boundary tracer values

            C0 = zeros(model_dims) # zero interior tracer to identify boundary source
            Jb_local = local_boundary_flux( f, C0, Fb)
            Jb = boundary_flux( f, C0, Fb)

            # check: filled with zeroes away from boundary?
            @test isequal(sum(Jb), sum(Jb_local))

            @testset "construct transport matrix" begin 
                # given Cb and mass circulation, solve for dC/dt
                Crand = rand(model_dims)

                # boundary flux is already expressed as a convergence        
                deldotJ = convergence(
                    advective_diffusive_flux(Crand, Fv))
                + boundary_flux(f, Crand, Fb)

                # ease the programming with a top-level driver function
                dCdt = tracer_tendency(Crand, f, Fv, Fb, Vol)
        
                # find A matrix.
                # If f = 0, q = 0, then dC/dt  = Ac
                A =  linear_probe(tracer_tendency, Crand, f, Fv, Fb, Vol)

                # view matrix in usual mathematical form
                Matrix(A)

                # probe for B (boundary matrix)
                dCdt_boundary = tracer_tendency(f, Crand, Fb, Vol)
                B =  linear_probe(tracer_tendency, f, Crand, Fb, Vol)
                Matrix(B)

                mass(Vol)

                # Find eigenvalues of A. 
                # destructuring via iteration
                μ, V = eigen(A)

                Tmax = maximum_timescale(μ)

                # water-mass fractions


                tmp = μ_matrix \ (V \ ustrip.(B))
                tmp = μ_matrix \ (V \ B)

                a = - real.(V * tmp)

                a = watermass_fraction(μ, V, B)
                Matrix(a)
                @test all(isapprox.(1.0,sum(a)))                
            end
        end
        
   end
    
end

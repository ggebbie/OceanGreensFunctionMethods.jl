using Revise
using OceanGreensFunctionMethods
using Distributions
using DimensionalData
using DimensionalData: @dim
using Unitful
using MultipliableDimArrays
using LinearAlgebra
# using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomadinSets, Unitful
using Test

#include(srcdir("config_units.jl")) # not sure why it doesn't work
include("../src/config_units.jl")

@testset "OceanGreensFunctionMethods.jl" begin

    @testset "tracer inverse Gaussian distribution" begin
        # compare a Tracer Inverse Gaussian distribution
        # with geophysical arguments
        # to an Inverse Gaussian with statistical arguments
        Γ = 20.0 # mean
        Δ = 20.0 # width

        G1 = TracerInverseGaussian(Γ, Δ)

        # the Inverse Gaussian distribution in typical statistical toolboxes
        # has different input arguments
        μ = Γ # the same, but different notation
        λ = shape(G1) # the shape parameter is the second argument 
        
        G2 = InverseGaussian(μ, λ)

        # both Inverse Gaussians should be the same
        @test isequal(shape(G1),shape(G2))
        @test !isequal(params(G1),params(G2))
        @test isequal(partype(G1), partype(G2)) # parameter type
        @test isequal(mean(G1), mean(G2))
        @test isequal(var(G1), var(G2))
    end

    @testset "pedagogical tracer box model" begin

        # define grid
        model_dims = model_dimensions()
        Ny, Nz = size(model_dims) # size in each dimension
        Nb = Ny * Nz # number of boxes
        
        #Vol_uniform = 1e16m^3 |> km^3 # uniform value of volume for all boxes
        Vol_uniform = 300.0Sv*yr |> km^3 # uniform value of volume for all boxes
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
            boundary_dims = boundary_dimensions()
           
            Fb = DimArray(hcat([20Sv, 10Sv]), boundary_dims) # boundary flux
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

                @testset "matrix exponential" begin
                    dt = 0.1yr
                    Matrix(μ)
                    Matrix(μ*dt)
                    exp(Matrix(μ*dt))
                    Matrix(A)
                    matexp = MultipliableDimArray( exp(Matrix(μ*dt)), dims(μ), dims(μ))
                    t1 =  real.( V * (matexp * (V\C))) # matlab code has right divide (?)

                    # mostly handled by MultipliableDimArrays 
                    #eAt = MultipliableDimArray(exp(Matrix(A*dt)),dims(A),dims(A))
                    eAt = exp(A*dt)
                    t2 = real.( eAt*C) # matlab code has right divide (?)
                    t3 = vec(t1) - vec(t2)
                    @test maximum(abs.(t3)) < 1e-8
                end
                Tmax = maximum_timescale(μ)

                # water-mass fractions
                a = watermass_fraction(μ, V, B, alg=:forward)
                Matrix(a)
                @test all(isapprox.(1.0,sum(a)))                

                a_adjoint = watermass_fraction(μ, V, B, alg=:adjoint)
                Matrix(a_adjoint)
                @test all(isapprox.(1.0,sum(a)))                

                a_residence = watermass_fraction(μ, V, B, alg=:residence)
                Matrix(a_residence)
                @test all(isapprox.(1.0,sum(a)))                

                Γ = mean_age(μ, V, B, alg=:forward)
                @test all(Γ .≥ 0.0yr)

                Γ_adjoint = mean_age(μ, V, B, alg=:adjoint)
                @test all(Γ_adjoint .≥ 0.0yr)

                # very similar values; is this correct?
                Δ = ttd_width(μ, V, B)
                @test 90yr < Δ[2,2] < 91yr # compare to MATLAB point value
                @test all(Δ .≥ 0.0yr)

                @testset "green's function" begin
                    Δτ = 0.25yr
                    τ = 0yr:Δτ:2000yr
                    ttest = 1.0yr
                    G(t) = greens_function(t,A) # a closure that captures A
                    @test all(Matrix(G(ttest)) .≥ 0.0)

                    # add test: normalization of Green's function
                    
                    G′(t) = boundary_propagator(t,A,B, alg=:forward)
                    @test all(Matrix(G′(ttest)) .≥ 0.0/yr)

                    # † is invalid in Julia as an identifier 
                    G′dagger(t) = boundary_propagator(t,A,B, alg=:adjoint)
                    @test all(Matrix(G′dagger(ttest)) .≥ 0.0/yr)

                    # residence times
                    # numerical values quite different from MATLAB
                    a_RTD = residence_time_watermass_fractions(μ, V, B)
                    @test isapprox(sum(Matrix(a_RTD)[:]),1.0) # fails due to missing complex conjugate
                end

                @testset "read tracer histories" begin

                    BD = read_tracer_histories()
                    tracername = :CFC11NH
                    box2_box1_ratio = 0.75
                    source_history_func(x) =  tracer_source_history(x, tracername, BD, box2_box1_ratio)
                
                    tt = 1973.0yr
                    source_history_func(tt)

                    ti = 1980.0yr
                    tf = 1981.0yr
                    source_history_func(tf)
                    tester = integrate_forcing(ti, tf, μ, V, B, source_history_func)

                    # goal: source_history(t,tracerHistory,radio_tracer,Tracer.(radio_tracer).box2_box1_ratio) ;

                    C₀ = zeros(model_dims)
                    tlist = (1980.0:1981.0)yr
                    # tmp = Array{DimArray}(undef,size(tlist))
                    # Cevolve = DimArray(tmp,Ti(tlist))

                    Cevolve = evolve_concentration(C₀, A, B, tlist, source_history_func; halflife = nothing)

                    sss =  [Cevolve[t][3,1] for t in eachindex(tlist)]
                end
            end
        end
    end
end


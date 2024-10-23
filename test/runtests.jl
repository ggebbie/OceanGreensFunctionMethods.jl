using Revise
using OceanGreensFunctionMethods
using AlgebraicArrays
using Distributions
using DimensionalData
using DimensionalData: @dim
using Unitful
using LinearAlgebra
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
        #Vol = DimArray(fill(Vol_uniform, Ny, Nz), model_dims)
        #Vol = VectorArray(fill(Vol_uniform, model_dims))
        Vol = fill(Vol_uniform, model_dims, :VectorArray)
 
        Ψ_abyssal = 20Sv
        Fv_abyssal = abyssal_overturning(Ψ_abyssal, model_dims) # volume fluxes

        Ψ_intermediate = 10Sv
        Fv_intermediate = intermediate_overturning(Ψ_intermediate, model_dims) # volume fluxes

        Fv_exchange = 5Sv
        Fv_diffusion = vertical_diffusion(Fv_exchange, model_dims) # volume fluxes

        Fv = Fv_abyssal + Fv_intermediate + Fv_diffusion
        C = ones(model_dims, :VectorArray)

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

        # boundary exchange: define the locations affected by boundary fluxes
        boundary_dims = boundary_dimensions()

        # AlgebraicArray takes an object with dimensional info and makes it a VectorArray or MatrixArray
        Fb = AlgebraicArray([20Sv, 10Sv], boundary_dims) # boundary flux
        f = ones(boundary_dims, :VectorArray) # boundary tracer values

        C0 = zeros(model_dims, :VectorArray) # zero interior tracer to identify boundary source
        Jb_local = local_boundary_flux( f, C0, Fb)
        Jb = boundary_flux( f, C0, Fb)

        # check: filled with zeroes away from boundary?
        @test isequal(sum(Jb), sum(Jb_local))

        @testset "construct transport matrix" begin 
            # given Cb and mass circulation, solve for dC/dt
            # boundary flux is already expressed as a convergence        
            deldotJ = convergence(
                advective_diffusive_flux(C0, Fv))
            + boundary_flux(f, C0, Fb)

            # ease the programming with a top-level driver functions
            dCdt = tracer_tendency(C0, f, Fv, Fb, Vol)
            dCdt_boundary = tracer_tendency(f, C0, Fb, Vol)
            dCdt_radioactive = tracer_tendency(C0, 269yr)

            # should be true, but is not
            typeof(tracer_tendency(C0, f, Fv, Fb, Vol) - tracer_tendency(C0, f, Fv, Fb, Vol)) ==
            typeof(tracer_tendency(C0, f, Fv, Fb, Vol))

            # find A matrix.
            # If f = 0, q = 0, then dC/dt  = Ac

            A =  linear_probe(tracer_tendency, C0, f, Fv, Fb, Vol)

            # probe for B (boundary matrix)
            B =  linear_probe(tracer_tendency, f, C0, Fb, Vol)

            # probe for radioactive adjustment to matrix
            Aλ =  linear_probe(tracer_tendency, C0, 269yr)

            # Find eigenvalues of A. 
            # destructuring via iteration
            μ, V = eigen(A)

            Tmax = maximum_timescale(μ)

            @testset "matrix exponential" begin
                dt = 0.1yr
                eAt = exp(A*dt)
                t2 = real.( eAt*C) # matlab code has right divide (?)
                @test maximum(t2) ≤ 1.0
            end

            @testset "global TTD" begin
                𝒢(t) = global_ttd(t,A,B) # type \scr + G + TAB
                ttd_global = 𝒢(1yr)[At("High latitudes"),At("Thermocline")]
                @test ttd_global ≥ 0.0/yr
            end
            
            @testset "water masses" begin

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

            end


            @testset "mean and ideal ages" begin
            
                Γ = mean_age(μ, V, B, alg=:forward)
                @test all(Γ .≥ 0.0yr)

                Γ_adjoint = mean_age(μ, V, B, alg=:adjoint)
                @test all(Γ_adjoint .≥ 0.0yr)

                Γ_residence = mean_age(μ, V, B, alg=:residence)
                @test 258yr < Γ_residence < 259yr

                Γ_ideal = ideal_age(A, B, alg= :forward)
                # should be identical to mean age
                                    
                Γ_ideal_adjoint = ideal_age(A, B, alg= :adjoint)

                # test one random location (could test all locations)
                itest = rand(1:Ny)
                jtest = rand(1:Nz)
                @test isapprox(Γ_ideal[itest,jtest], Γ[itest,jtest], rtol = 1e-3)
                @test isapprox(Γ_ideal_adjoint[itest,jtest], Γ_adjoint[itest,jtest], rtol = 1e-3)
            end

            @testset "TTD width" begin

                Δ = ttd_width(μ, V, B)
                @test 90yr < Δ[2,2] < 91yr # compare to MATLAB point value
                @test all(Δ .≥ 0.0yr)

                Δ_adjoint = ttd_width(μ, V, B, alg=:adjoint)
                @test 90yr < Δ_adjoint[2,2] < 91yr # compare to MATLAB point value
                @test all(Δ_adjoint .≥ 0.0yr)

                Δ_residence = ttd_width(μ, V, B, alg=:residence)
                @test 129yr < Δ_residence < 130yr # compare to MATLAB point value

            end

            @testset "green's function" begin
                Δτ = 0.25yr
                τ = 0yr:Δτ:2000yr
                ttest = 1.0yr
                G(t) = greens_function(t, A) # a closure that captures A
                @test all(Matrix(G(ttest)) .≥ 0.0)

                # missing test: normalization of Green's function
                    
                G′(t) = boundary_propagator(t, A, B, alg=:forward)
                @test all(Matrix(G′(ttest)) .≥ 0.0/yr)

                # † is invalid in Julia as an identifier 
                G′dagger(t) = boundary_propagator(t, A, B, alg=:adjoint)
                @test all(Matrix(G′dagger(ttest)) .≥ 0.0/yr)

                𝒢(t) = global_ttd(t, A, B, alg=:forward)

                𝒢dagger(t) = global_ttd(t, A, B, alg=:adjoint)
                𝒢dagger(1yr)

                RTD(t) = residence_time(t,A,B)
                RTD(1yr)
                    
                # residence times
                # check: numerical values different from MATLAB?
                a_residence = watermass_fraction(μ, V, B, alg=:residence)
                @test isapprox(sum(Matrix(a_residence)),1.0) 
            end

            @testset "path density" begin
                Φ(τ) = OceanGreensFunctionMethods.phi_function(τ, μ) # a useful closure
                Matrix(Φ(10yr))
                # missing test for proper normalization

                mbox = "High latitudes"
                vbox = "Thermocline"
                D_mat = AlgebraicArray(zeros(length(V), length(V)),model_dimensions(),model_dimensions())
                D_mat[At(mbox),At(vbox)][At(mbox),At(vbox)] = 1 
                D_mat_overline = V \ D_mat * V

                # check for element-by-element product to simplify 
                elemental_product = OceanGreensFunctionMethods.hadamard(D_mat_overline,Φ(1yr))
                pd = path_density(μ, V, B, 1yr, mbox, vbox)
                @test all( Matrix(pd) .≥ 0.0/yr)
            end

            @testset "read tracer histories" begin

                BD = read_transient_tracer_histories()
                tracername = :CFC11NH
                box2_box1_ratio = 0.75

                @test tracer_source_history(1990yr,
                    tracername,
                    box2_box1_ratio,
                    BD) isa VectorDimArray
                    
                source_history_func1(t) =  tracer_source_history(t,
                    tracername,
                    box2_box1_ratio,
                    BD,
                )
                    
                tt = 1973.0yr
                source_history_func1(tt)

                ti = 1980.0yr
                tf = 1981.0yr
                source_history_func1(tf)
                func_test(t) = OceanGreensFunctionMethods.forcing_integrand(t, tf, μ, V, B, source_history_func1)
                tester = integrate_forcing(ti, tf, μ, V, B, source_history_func1) # does it run?

                C₀ = zeros(model_dims, :VectorArray)
                tlist = (1980.0:1981.0)yr
                Cevolve = evolve_concentration(C₀, A, B, tlist, source_history_func1; halflife = nothing)
                Ct =  [Cevolve[t][3,1] for t in eachindex(tlist)]
                @test Ct[end] > Ct[begin] 

                # argon-39
                tracername = :argon39
                box2_box1_ratio = 1 
                source_history_func2(t) =  tracer_source_history(t,
                    tracername,
                    box2_box1_ratio,
                )
                tt = 1973.0yr
                # always returns 1 
                @test isequal(first(source_history_func2(2000yr*randn())),1.0)

                # iodine-129
                BD_iodine129 = read_iodine129_history()
                tracername = :iodine129
                box2_box1_ratio = 0.25

                tracer_source_history(1990yr,
                    tracername,
                    box2_box1_ratio,
                    BD_iodine129)
                    
                source_history_func3(t) =  tracer_source_history(t,
                    tracername,
                    box2_box1_ratio,
                    BD_iodine129,
                )
                    
                tt = 1873.0yr
                @test source_history_func3(tt) isa VectorDimArray
                    
            end

        end
    end
end


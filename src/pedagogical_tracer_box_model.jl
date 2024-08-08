# follows the design of Tom Haine's pedagaogical box model

kg = u"kg"
m  = u"m"
yr = u"yr"
Tg = u"Tg"
s  = u"s"

struct Fluxes{T,N} 
    poleward::DimArray{T,N}
    equatorward::DimArray{T,N}
    up::DimArray{T,N}
    down::DimArray{T,N}
end

+(F1::Fluxes, F2::Fluxes) = Fluxes(
    F1.poleward .+ F2.poleward,
    F1.equatorward .+ F2.equatorward,
    F1.up .+ F2.up,
    F1.down .+ F2.down)

dims(F::Fluxes) = dims(F.poleward)

meridional_names() = ["1 High latitudes", "2 Mid-latitudes", "3 Low latitudes"]
vertical_names() = ["1 Thermocline", "2 Deep", "3 Abyssal"]

function abyssal_overturning(Ψ,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Ψ)
    Fv_poleward = zeros(model_dims)*Fv_units
    Fv_equatorward = zeros(model_dims)*Fv_units
    Fv_up = zeros(model_dims)*Fv_units
    Fv_down = zeros(model_dims)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_poleward[At("3 Low latitudes"),At("1 Thermocline")] = Ψ 
    Fv_poleward[At("2 Mid-latitudes"),At("1 Thermocline")] = Ψ 

    Fv_equatorward[At("2 Mid-latitudes"),At("3 Abyssal")] = Ψ 
    Fv_equatorward[At("1 High latitudes"),At("3 Abyssal")] = Ψ 

    Fv_up[At("3 Low latitudes"),At("3 Abyssal")] = Ψ 
    Fv_up[At("3 Low latitudes"),At("2 Deep")] = Ψ 

    Fv_down[At("1 High latitudes"),At("1 Thermocline")] = Ψ 
    Fv_down[At("1 High latitudes"),At("2 Deep")] = Ψ 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

function intermediate_overturning(Ψ,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Ψ)
    Fv_poleward = zeros(model_dims)*Fv_units
    Fv_equatorward = zeros(model_dims)*Fv_units
    Fv_up = zeros(model_dims)*Fv_units
    Fv_down = zeros(model_dims)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_poleward[At("2 Mid-latitudes"),At("3 Abyssal")] = Ψ 
    Fv_poleward[At("3 Low latitudes"),At("3 Abyssal")] = Ψ 

    Fv_equatorward[At("1 High latitudes"),At("2 Deep")] = Ψ 
    Fv_equatorward[At("2 Mid-latitudes"),At("2 Deep")] = Ψ 

    Fv_up[At("1 High latitudes"),At("3 Abyssal")] = Ψ 
    Fv_down[At("3 Low latitudes"),At("2 Deep")] = Ψ 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

function vertical_diffusion(Fv_exchange,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Fv_exchange)
    Fv_poleward = zeros(model_dims)*Fv_units
    Fv_equatorward = zeros(model_dims)*Fv_units
    Fv_up = zeros(model_dims)*Fv_units
    Fv_down = zeros(model_dims)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_up[:,At(["3 Abyssal","2 Deep"])] .= Fv_exchange 
    Fv_down[:,At(["1 Thermocline","2 Deep"])] .= Fv_exchange 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

advective_diffusive_flux(C::DimArray, Fv::DimArray ; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

advective_diffusive_flux(C::DimArray, Fv::Fluxes ; ρ = 1035kg/m^3) =
    Fluxes(
        advective_diffusive_flux(C, Fv.poleward, ρ=ρ),
        advective_diffusive_flux(C, Fv.equatorward, ρ=ρ),
        advective_diffusive_flux(C, Fv.up, ρ=ρ),
        advective_diffusive_flux(C, Fv.down, ρ=ρ)
    )

mass(V ; ρ = 1035kg/m^3) = ρ * V .|> u"Zg"

function convergence(J::Fluxes)

    # all the fluxes leaving a box
    deldotJ = -( J.poleward + J.equatorward + J.up + J.down)

    #poleward flux entering
    deldotJ[At(["2 Mid-latitudes","1 High latitudes"]),:] .+=
        J.poleward[At(["3 Low latitudes","2 Mid-latitudes"]),:]

    #equatorward flux entering
    deldotJ[At(["3 Low latitudes","2 Mid-latitudes"]),:] .+=
        J.equatorward[At(["2 Mid-latitudes","1 High latitudes"]),:]

    # upward flux entering
    deldotJ[:,At(["1 Thermocline","2 Deep"])] .+=
        J.up[:,At(["2 Deep","3 Abyssal"])]

    # downward flux entering
    deldotJ[:,At(["2 Deep","3 Abyssal"])] .+=
        J.down[:,At(["1 Thermocline","2 Deep"])]

    return deldotJ 
end

mass_convergence(Fv) = convergence(advective_diffusive_flux(ones(dims(Fv)), Fv))

function local_boundary_flux(f::DimArray, C::DimArray, Fb::DimArray)
    ΔC = f - C[DimSelectors(f)] # relevant interior tracer difference from boundary value
    return Jb = advective_diffusive_flux(ΔC, Fb)
end

function boundary_flux(f::DimArray, C::DimArray, Fb::DimArray)
    Jlocal = local_boundary_flux(f, C, Fb)
    Jb = unit(first(Jlocal)) * zeros(dims(C)) # pre-allocate
    Jb[DimSelectors(f)] += Jlocal # transfer J at boundary locations onto global grid
    return Jb
end

# tracer_tendency(C::DimArray{T,N}, Fv::Fluxes{T,N}) where T <: Number where N = 
#         convergence(advective_diffusive_flux(C, Fv))

tracer_tendency(
    C::DimArray{<:Number,N},
    f::DimArray{<:Number,N},
    Fv::Fluxes{<:Number,N},
    Fb::DimArray{<:Number,N},
    V::DimArray{<:Number,N}) where N = 
    ((convergence(advective_diffusive_flux(C, Fv)) +
    boundary_flux(f, C, Fb)) ./
    mass(V)) .|> yr^-1 

    # for use with finding B boundary matrix
tracer_tendency(
    f::DimArray{<:Number,N},
    C::DimArray{<:Number,N},
    Fb::DimArray{<:Number,N},
    V::DimArray{<:Number,N}) where N = 
    (boundary_flux(f, C, Fb) ./
    mass(V)) .|> yr^-1 

"""
function linear_probe(x₀,M)

    Probe a function to determine its linear response in matrix form.
    Assumes units are needed and available.
    A simpler function to handle cases without units would be nice.

    funk:: function to be probed
    x:: input variable
    args:: the arguments that follow x in `funk`
"""
function linear_probe(funk::Function,C::DimArray{T,N},args...) where T <: Number where N
    #function convolve(P::DimArray{T},M) where T<: AbstractDimArray

    dCdt0 = funk(C, args...)
    Trow = typeof(dCdt0)
    A = Array{Trow}(undef,size(C))

    for i in eachindex(C)
        C[i] += 1.0*unit(first(C))
        # remove baseline if not zero
        # not strictly necessary for linear system 
        A[i] = funk(C, args...) - dCdt0
        C[i] -= 1.0*unit(first(C)) # necessary?
    end
    return DimArray(A, dims(C))
end

function eigen(A::DimArray{<:DimArray})

    uniform(A) ? uA = unit(first(first(A))) : error("No eigendecomposition for a non-uniform matrix")
    A_matrix = MultipliableDimArrays.Matrix(A)
    F = eigen(ustrip.(A_matrix))
    values = uA * F.values

    eigen_dims = Eigenmode(1:size(A_matrix,2))
    model_dims = dims(A)
    vectors = MultipliableDimArray(F.vectors,
            model_dims, eigen_dims)    

    println("return vals")
    return values, vectors
    # ideally, would return an Eigen factorization, in spirit like:
    #    return Eigen(QuantityArray(F.values, dimension(A)), F.vectors)
end

allequal(x) = all(y -> y == first(x), x)

# eigenstructure only exists if A is uniform
function uniform(A::DimArray{<:DimArray})
    ulist = unit.(MultipliableDimArrays.Matrix(A))
    return allequal(ulist)
end

maximum_timescale(μ) = -1/real.(μ)[end]
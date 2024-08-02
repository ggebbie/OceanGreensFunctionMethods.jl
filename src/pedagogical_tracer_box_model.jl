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

tracer_flux(Fv::DimArray, C::DimArray; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

tracer_flux(Fv::Fluxes, C::DimArray; ρ = 1035kg/m^3) =
    Fluxes(
        tracer_flux( Fv.poleward, C, ρ=ρ),
        tracer_flux( Fv.equatorward, C, ρ=ρ),
        tracer_flux( Fv.up, C, ρ=ρ),
        tracer_flux( Fv.down, C, ρ=ρ)
    )

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

mass_convergence(Fv) = convergence(tracer_flux(Fv,ones(dims(Fv))))

#boundary_flux(Fv::DimArray, C::DimArray; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

function boundary_flux(Fb::DimArray,Cb::DimArray,C::DimArray)
    ΔC = Cb - C[DimSelectors(Cb)] # relevant interior tracer difference from boundary value
    Jb = tracer_flux(Fb, ΔC)
end

function apply_boundary_flux(Fb::DimArray,Cb::DimArray,C::DimArray)
    Jb = 0.0 * similar(C)
    Jb[DimSelectors(Cb)] += boundary_flux(Fb,Cb,C)
end

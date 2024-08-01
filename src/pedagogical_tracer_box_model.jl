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

function abyssal_overturning(Ψ,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Ψ)
    Fv_poleward = zeros(model_dims)*Fv_units
    Fv_equatorward = zeros(model_dims)*Fv_units
    Fv_up = zeros(model_dims)*Fv_units
    Fv_down = zeros(model_dims)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_poleward[At("Low latitudes"),At(:Thermocline)] = Ψ 
    Fv_poleward[At("Mid-latitudes"),At(:Thermocline)] = Ψ 

    Fv_equatorward[At("Mid-latitudes"),At(:Abyssal)] = Ψ 
    Fv_equatorward[At("High latitudes"),At(:Abyssal)] = Ψ 

    Fv_up[At("Low latitudes"),At(:Abyssal)] = Ψ 
    Fv_up[At("Low latitudes"),At(:Deep)] = Ψ 

    Fv_down[At("High latitudes"),At(:Thermocline)] = Ψ 
    Fv_down[At("High latitudes"),At(:Deep)] = Ψ 

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
    Fv_poleward[At("Mid-latitudes"),At(:Abyssal)] = Ψ 
    Fv_poleward[At("Low latitudes"),At(:Abyssal)] = Ψ 

    Fv_equatorward[At("High latitudes"),At(:Deep)] = Ψ 
    Fv_equatorward[At("Mid-latitudes"),At(:Deep)] = Ψ 

    Fv_up[At("High latitudes"),At(:Abyssal)] = Ψ 
    Fv_down[At("Low latitudes"),At(:Deep)] = Ψ 

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
    Fv_up[:,At([:Abyssal,:Deep])] .= Fv_exchange 
    Fv_down[:,At([:Thermocline,:Deep])] .= Fv_exchange 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

tracer_flux(Fv::DimArray, C::DimArray; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

function tracer_flux(Fv::Fluxes, C::DimArray; ρ = 1035kg/m^3)
return Fluxes(
    tracer_flux( Fv.poleward, C, ρ=ρ),
    tracer_flux( Fv.equatorward, C, ρ=ρ),
    tracer_flux( Fv.up, C, ρ=ρ),
    tracer_flux( Fv.down, C, ρ=ρ))

end

function tracer_flux_convergence(J::Fluxes)

    # all the fluxes leaving a box
    deldotJ = -( J.poleward + J.equatorward + J.up + J.down)

    #poleward flux entering
    deldotJ[At(["Mid-latitudes","High latitudes"]),:] .+=
        J.poleward[At(["Low latitudes","Mid-latitudes"]),:]

    #equatorward flux entering
    deldotJ[At(["Low latitudes","Mid-latitudes"]),:] .+=
        J.equatorward[At(["Mid-latitudes","High latitudes"]),:]

    # upward flux entering
    deldotJ[:,At([:Thermocline,:Deep])] .+=
        J.up[:,At([:Deep,:Abyssal])]

    # downward flux entering
    deldotJ[:,At([:Deep,:Abyssal])] .+=
        J.down[:,At([:Thermocline,:Deep])]

    return deldotJ 
end

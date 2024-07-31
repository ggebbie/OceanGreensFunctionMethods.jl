# follows the design of Tom Haine's pedagaogical box model

kg = u"kg"
m  = u"m"
yr = u"yr"
Tg = u"Tg"
s  = u"s"

struct VolumeFlux{T,N} where T <: Number where N 
    poleward::DimArray{T,N}
    equatorward::DimArray{T,N}
    up::DimArray{T,N}
    down::DimArray{T,N}
    dims::Tuple
end

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

    return Fv_poleward, Fv_equatorward, Fv_up, Fv_down
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
    Fv_poleward[At("High latitudes"),At(:Abyssal)] = Ψ 

    Fv_equatorward[At("Low latitudes"),At(:Deep)] = Ψ 
    Fv_equatorward[At("Mid-latitudes"),At(:Deep)] = Ψ 

    Fv_up[At("High latitudes"),At(:Abyssal)] = Ψ 
    Fv_down[At("Low latitudes"),At(:Deep)] = Ψ 

    return Fv_poleward, Fv_equatorward, Fv_up, Fv_down
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

    return Fv_poleward, Fv_equatorward, Fv_up, Fv_down
end

tracer_flux(Fv, C; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

#end

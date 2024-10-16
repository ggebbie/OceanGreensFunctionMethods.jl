# follows the design of Tom Haine's pedagaogical box model

# define units
kg = u"kg"
m  = u"m"
yr = u"yr"
Tg = u"Tg"
s  = u"s"
pmol = u"pmol"
fmol = u"fmol"
nmol = u"nmol"

# define dimensional labels
@dim Tracer "tracer"
@dim Meridional "meridional location"
@dim Vertical "vertical location"
@dim Global "global quantity"

# define a structure for 2D fluxes in a yz domain
struct Fluxes{T,A<:AbstractArray{T}} 
    poleward::A
    equatorward::A
    up::A
    down::A
end

+(F1::Fluxes, F2::Fluxes) = Fluxes(
    F1.poleward + F2.poleward,
    F1.equatorward + F2.equatorward,
    F1.up + F2.up,
    F1.down + F2.down)

dims(F::Fluxes) = dims(F.poleward)

meridional_names() = ["High latitudes", "Mid-latitudes", "Low latitudes"]
vertical_names() = ["Thermocline", "Deep", "Abyssal"]

"""
    model_dimensions()

Define labels for the model's physical dimensions, as well as labels for the box names. Use the format of `DimensionalData.jl`. Permits numerical quantities to be bundled with their meta-data. Dimensions are `Unordered` to avoid issues related to the alphabetical order.
"""
model_dimensions() = (Meridional(meridional_names(); order=DimensionalData.Unordered()),
    Vertical(vertical_names(); order=DimensionalData.Unordered())) 

"""
    boundary_dimensions()

Define labels for the boundary's physical dimensions, as well as labels for the box names, consistently with the model dimensions. Use the format of `DimensionalData.jl`. Permits numerical quantities to be bundled with their meta-data. Dimensions are `Unordered` to avoid issues related to the alphabetical order.
"""
boundary_dimensions() = (Meridional(meridional_names()[1:2]; order=DimensionalData.Unordered()),
    Vertical([vertical_names()[1]]; order=DimensionalData.Unordered())) 

# function Base.zeros(model_dims, type::Symbol)
#     if type == :Fluxes
#         Fv_zeros = zeros(model_dims, :VectorArray)*Fv_units
#         return Fluxes(Fv_zeros, Fv_zeros, Fv_zeros, Fv_zeros)
#     else
#         error("zeros not implemented for this type")
#     end
# end
    
"""
    abyssal_overturning(Ψ,model_dims)

Set volume flux, Ψ, in an abyssal overturning loop that satisfies the conservation of volume. Return a structure of `Fluxes`.
"""
function abyssal_overturning(Ψ,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Ψ)
    Fv_poleward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_equatorward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_up = zeros(model_dims, :VectorArray)*Fv_units
    Fv_down = zeros(model_dims, :VectorArray)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_poleward[At("Low latitudes"),At("Thermocline")] = Ψ 
    Fv_poleward[At("Mid-latitudes"),At("Thermocline")] = Ψ 

    Fv_equatorward[At("Mid-latitudes"),At("Abyssal")] = Ψ 
    Fv_equatorward[At("High latitudes"),At("Abyssal")] = Ψ 

    Fv_up[At("Low latitudes"),At("Abyssal")] = Ψ 
    Fv_up[At("Low latitudes"),At("Deep")] = Ψ 

    Fv_down[At("High latitudes"),At("Thermocline")] = Ψ 
    Fv_down[At("High latitudes"),At("Deep")] = Ψ 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

"""
    intermediate_overturning(Ψ,model_dims)

Set the volume flux, Ψ, in an intermediate overturning loop that satisfies the conservation of volume. Return a structure of `Fluxes`.
"""
function intermediate_overturning(Ψ,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Ψ)
    Fv_poleward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_equatorward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_up = zeros(model_dims, :VectorArray)*Fv_units
    Fv_down = zeros(model_dims, :VectorArray)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    Fv_poleward[At("Mid-latitudes"),At("Abyssal")] = Ψ 
    Fv_poleward[At("Low latitudes"),At("Abyssal")] = Ψ 

    Fv_equatorward[At("High latitudes"),At("Deep")] = Ψ 
    Fv_equatorward[At("Mid-latitudes"),At("Deep")] = Ψ 

    Fv_up[At("High latitudes"),At("Abyssal")] = Ψ 
    Fv_down[At("Low latitudes"),At("Deep")] = Ψ 

    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end


"""
    vertical_diffusion(Fv_exchange,model_dims)

Set vertical diffusive-like exchange flux `Fv_exchange`. Return a structure of `Fluxes`.
"""
function vertical_diffusion(Fv_exchange,model_dims)

    # pre-allocate volume fluxes with zeros with the right units
    Fv_units = unit(Fv_exchange)
    Fv_poleward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_equatorward = zeros(model_dims, :VectorArray)*Fv_units
    Fv_up = zeros(model_dims, :VectorArray)*Fv_units
    Fv_down = zeros(model_dims, :VectorArray)*Fv_units

    # set fluxes manually
    # fluxes organized according to (upwind) source of flux
    # missing proper broadcast for VectorDimArray: add `parent` below
    Fv_up[:,At(["Abyssal","Deep"])] .= Fv_exchange 
    Fv_down[:,At(["Thermocline","Deep"])] .= Fv_exchange 
    #parent(Fv_up)[:,At(["Abyssal","Deep"])] .= Fv_exchange 
    #parent(Fv_down)[:,At(["Thermocline","Deep"])] .= Fv_exchange 
    
    return Fluxes(Fv_poleward, Fv_equatorward, Fv_up, Fv_down)
end

"""
    advective_diffusive_flux(C, Fv; ρ)

Advective-diffusive flux of tracer `C` given volume fluxes `Fv` and optional density `ρ`.

# Arguments
- `C::DimArray`: tracer distribution
- `Fv::DimArray`: volume fluxes
- `ρ::Number=1035kg/m^3`: uniform density
# Returns
- `Fc::DimArray`: tracer flux
"""
advective_diffusive_flux(C::VectorDimArray, Fv::VectorDimArray ; ρ = 1035kg/m^3) = ρ * (Fv .* C) .|> Tg/s

"""
    advective_diffusive_flux(C, Fv; ρ)

Advective-diffusive flux of tracer `C` given volume fluxes `Fv` and optional density `ρ`.

# Arguments
- `C::DimArray`: tracer distribution
- `Fv::Fluxes`: volume fluxes
- `ρ::Number=1035kg/m^3`: uniform density
# Returns
- `Fc::Fluxes`: tracer flux
"""
advective_diffusive_flux(C::VectorDimArray, Fv::Fluxes ; ρ = 1035kg/m^3) =
    Fluxes(
        advective_diffusive_flux(C, Fv.poleward, ρ=ρ),
        advective_diffusive_flux(C, Fv.equatorward, ρ=ρ),
        advective_diffusive_flux(C, Fv.up, ρ=ρ),
        advective_diffusive_flux(C, Fv.down, ρ=ρ)
    )

"""    
    mass(V; ρ)

Seawater mass derived from the volume `V` and an optional input of density `ρ`.
"""
mass(V; ρ = 1035kg/m^3) = ρ * V .|> u"Zg"

"""
    convergence(J)

Convergence of fluxes `J` of type `Fluxes`.
This is a computational methods that depends on proper slices and broadcasting
and thus currently requires using `parent` on the left hand side below.
"""
function convergence(J::Fluxes{T,A}) where {T, A <: VectorDimArray{T}}

    # all the fluxes leaving a box
    deldotJ = -( J.poleward + J.equatorward + J.up + J.down)

    #poleward flux entering
    parent(deldotJ)[At(["Mid-latitudes","High latitudes"]),:] .+=
       J.poleward[At(["Low latitudes","Mid-latitudes"]),:]

    #equatorward flux entering
    parent(deldotJ)[At(["Low latitudes","Mid-latitudes"]),:] .+=
        J.equatorward[At(["Mid-latitudes","High latitudes"]),:]

    # upward flux entering
    parent(deldotJ)[:,At(["Thermocline","Deep"])] .+=
        J.up[:,At(["Deep","Abyssal"])]

    # downward flux entering
    parent(deldotJ)[:,At(["Deep","Abyssal"])] .+=
        J.down[:,At(["Thermocline","Deep"])]

    return deldotJ 
end

"""
    mass_convergence(Fv) 

Convergence of volume derived from a field of volume fluxes `Fv`, translated into a mass flux convergence with the assumption of uniform density.  
"""
mass_convergence(Fv) = convergence(advective_diffusive_flux( ones(dims(Fv), :VectorArray), Fv))

function local_boundary_flux(f::VectorDimArray,
    C::VectorDimArray,
    Fb::VectorDimArray)
    
    ΔC = f - C[DimSelectors(f)] # relevant interior tracer difference from boundary value
    return advective_diffusive_flux(ΔC, Fb)
end

"""
    boundary_flux(f::DimArray, C::DimArray, Fb::DimArray)

Convergence or net effect of boundary fluxes.

# Arguments
- `f::DimArray`: Dirichlet boundary condition
- `C::DimArray`: tracer distribution 
- `Fb::DimArray`: boundary exchange volume flux
# Returns
- `Jb::Fluxes`: boundary tracer flux
"""
function boundary_flux(f::VectorDimArray, C::VectorDimArray, Fb::VectorDimArray)
    Jlocal = local_boundary_flux(f, C, Fb)
    Jb = unit(first(Jlocal)) * VectorArray(zeros(dims(C))) # pre-allocate
    Jb[DimSelectors(f)] += Jlocal # transfer J at boundary locations onto global grid
    return Jb
end

"""
    radioactive_decay(C, halflife)

Radioactive decay rate of tracer `C` with half life of `halflife`.
"""
radioactive_decay(C::Union{VectorArray,DimArray}, halflife::Number) = -(log(2)/halflife)*C 

"""
    tracer_tendency(C, f, Fv, Fb, V)

Tracer tendency ∂C/∂t for a tracer `C`, especially useful for finding a tracer transport matrix. 

# Arguments
- `C::VectorDimArray`: tracer distribution
- `f::VectorDimArray`: Dirichlet boundary condition
- `Fv::Fluxes`: volume fluxes
- `Fb::VectorDimArray`: boundary flux convergence
- `V::VectorDimArray`: box volume
# Returns
- `dCdt::VectorDimArray`: tracer tendency
"""
tracer_tendency(
    C::VDA1,
    f::VDA2,
    Fv::Fluxes{T,VDA3},
    Fb::VDA4,
    V::VDA5) where {T, VDA1 <: VectorDimArray,
        VDA2 <: VectorDimArray,
        VDA3 <: VectorDimArray,
        VDA4 <: VectorDimArray,
        VDA5 <: VectorDimArray} =
    ((convergence(advective_diffusive_flux(C, Fv)) +
                  boundary_flux(f, C, Fb)) ./
                  mass(V)) .|> yr^-1 
# computationally fast version of tracer tendency
tracer_tendency(
    C::DimArray{T},
    f::DimArray{T},
    Fv::Fluxes{T},
    Fb::DimArray{T},
    V::DimArray{T}) where T <: Number = 
    ((convergence(advective_diffusive_flux(C, Fv)) +
                  boundary_flux(f, C, Fb)) ./
                  mass(V)) .|> yr^-1 

"""
    tracer_tendency(f, C, Fv, Fb, V)

Tracer tendency ∂C/∂t for a boundary flux `f`, for use with finding B boundary matrix.

# Arguments
- `f::DimArray`: Dirichlet boundary condition
- `C::DimArray`: tracer distribution
- `Fv::Fluxes`: volume fluxes
- `Fb::Fluxes`: volume fluxes
- `V::DimArray`: box volume
# Returns
- `dCdt::DimArray`: tracer tendency
"""
tracer_tendency(
    # f::DimArray{<:Number,N},
    # C::DimArray{<:Number,N},
    # Fb::DimArray{<:Number,N},
    # V::DimArray{<:Number,N}) where N = 
    f::VDA2,
    C::VDA1,
    Fb::VDA4,
    V::VDA5) where {T, VDA1 <: VectorDimArray,
        VDA2 <: VectorDimArray,
        VDA4 <: VectorDimArray,
        VDA5 <: VectorDimArray} =
    (boundary_flux(f, C, Fb) ./
    mass(V)) .|> yr^-1 

# for use with finding A perturbation with radioactive decay
"""
    tracer_tendency(C)

Tracer tendency ∂C/∂t for the radioactive decay of a tracer `C` with half life `halflife`, for use with finding the radioactive contribution to a tracer transport matrix.

# Arguments
- `C::DimArray`: tracer distribution
- `halflife::Number`: radioactive half life
# Returns
- `dCdt::DimArray`: tracer tendency
"""
tracer_tendency(C::Union{VectorDimArray{T,N},DimArray{T,N}},
    halflife::T) where {T,N} =
    radioactive_decay(C, halflife) .|> yr^-1 

"""
    linear_probe(funk, x, args...)

Probe a function to determine its linear response in matrix form. Assumes units are needed and available. A simpler function to handle cases without units would be nice.

# Arguments
- `funk`: function to be probed
- `x`: input (independent) variable
- `halflife::Number`: radioactive half life
- `args`: the arguments that follow `x` in `funk`
# Returns
- `A::DimArray{DimArray}`: labeled transport information used in matrix operations 
"""
function linear_probe(funk::Function,C::DimArray{T,N},args...) where T <: Number where N

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

function linear_probe(funk::Function, C::VectorDimArray, args...)
    return MatrixArray(linear_probe(funk,parent(C),parent.(args...)))
end

allequal(x) = all(y -> y == first(x), x)

"""
    location_transient_tracer_histories()

URL of tracer source history file.
"""
location_transient_tracer_histories() = "https://github.com/ThomasHaine/Pedagogical-Tracer-Box-Model/raw/main/MATLAB/tracer_gas_histories.mat"

"""
    location_transient_tracer_histories()

URL of iodine-129 source history file.
"""
location_iodine129_history() = 
    "https://raw.githubusercontent.com/ThomasHaine/Pedagogical-Tracer-Box-Model/main/MATLAB/From%20John%20Smith/Input%20Function%20129I%20Eastern%20Norwegian%20Sea.csv"

"""
    read_iodine129_history()
"""
function read_iodine129_history()
    url = location_iodine129_history()
    !isdir(datadir()) && mkpath(datadir())
    filename = datadir("iodine129_history.nc")

    # allow offline usage if data already downloaded
    !isfile(filename) && Downloads.download(url,filename)

    # use CSV to open
    ds = CSV.File(filename)# , DataFrame)

    source_iodine = vcat(0.0,0.0,parse.(Float64,ds.Column2[4:end]))
    t_iodine = vcat(0.0,1957.0,parse.(Float64,ds.var"Atlantic Water in Eastern Norwegian Sea entering Barents Sea and West Spitzbergen Current"[4:end]))

    tracerdim = Tracer([:iodine129])
    timedim = Ti((t_iodine)yr)

    BD = zeros(tracerdim,timedim)
    BD[Tracer=At(:iodine129)] = source_iodine # note that this does NOT introduce a variable ``varname`` into scope

    return BD
end 

"""
    read_transient_tracer_histories()

Read transient tracer source histories and save as a `DimArray`. 
"""
function read_transient_tracer_histories()

    # download tracer history input (make this lazy)
    url = location_transient_tracer_histories()
    !isdir(datadir()) && mkpath(datadir())
    filename = datadir("tracer_histories.mat")

    # allow offline usage if data already downloaded
    !isfile(filename) && Downloads.download(url,datadir("tracer_histories.mat"))

    file = matopen(filename)

    # all matlab variables except Year
    varnames = Symbol.(filter(x -> x ≠ "Year", collect(keys(file))))
    tracerdim = Tracer(varnames)
    timedim = Ti(vec(read(file, "Year"))yr)

    BD = zeros(tracerdim,timedim)
    
    for v in varnames
        BD[Tracer=At(v)] = read(file, string(v)) # note that this does NOT introduce a variable ``varname`` into scope
    end
    close(file)
    return BD
end

tracer_units() = Dict(
    :CFC11NH => NoUnits,
    :CFC11SH => NoUnits,
    :CFC12NH => NoUnits,
    :CFC12SH => NoUnits,
    :argon39 => NoUnits,
    :iodine129 => NoUnits,
    :SF6NH => NoUnits,
    :SF6SH => NoUnits,
    :N2ONH => nmol/kg,
    :N2OSH => nmol/kg
    )

    # for non-transient tracers
function tracer_point_source_history(tracername)

    if tracername == :argon39
        return x -> 1.0 * tracer_units()[tracername]
    elseif tracername == :iodine129
        error("not implemented yet")
    end
end

"""
    tracer_point_source_history(tracername, BD)

Return a function that yields transient tracer source history (such as CFCs) at given time.

# Arguments
- `tracername`: name of tracer in source history file
- `BD::DimArray`: Dirichlet boundary condition compendium for many tracers
"""
function tracer_point_source_history(tracername, BD)

    tracer_timeseries = BD[Tracer=At(tracername)] * tracer_units()[tracername]

    return linear_interpolation(
        first(DimensionalData.index(dims(tracer_timeseries))),
        tracer_timeseries)
end

"""
    tracer_source_history(t, tracername, box2_box1_ratio, BD = nothing)

Return source history values for all boundary points.

# Arguments
- `t`: time
- `tracername`: name of tracer in source history file
- `box2_box1_ratio`: ratio of boundary condition value in Mid-latitudes to High Latitudes
- `BD::DimArray=nothing`: Dirichlet boundary condition compendium (optional)
"""
function tracer_source_history(t, tracername, box2_box1_ratio, BD = nothing)

    if tracername == :argon39
        source_func = tracer_point_source_history(tracername)
    else
        source_func = tracer_point_source_history(tracername, BD)
    end

    box1 = source_func(t)
    box2 = box2_box1_ratio * box1
    
    # replace this section with a function call.
    boundary_dims = boundary_dimensions()
    return DimArray(hcat([box1,box2]),boundary_dims)
end

"""
    evolve_concentration(C₀, A, B, tlist, source_history; halflife = nothing)

Integrate forcing vector over time to compute the concentration history. Find propagator by analytical expression using eigen-methods.

# Arguments
- `C₀`: initial tracer concentration
- `A`: tracer transport information used in matrix calculations
- `B`: boundary condition information used in matrix calculations
- `tlist`: list of times to save tracer concentration
- `source_history::Function`: returns Dirichlet boundary condition at a given time
- `halflife=nothing`: radioactive half life (optional)
"""
function evolve_concentration(C₀, A, B, tlist, source_history; halflife = nothing)

    if isnothing(halflife)
        μ, V = eigen(A)
    else
        Aλ =  linear_probe(tracer_tendency, C₀, halflife)
        μ, V = eigen(A+Aλ)
    end 

    # initial condition contribution
    Ci = deepcopy(C₀)

    # forcing contribution
    Cf = zeros(dims(C₀))

    # total
    C = DimArray(Array{VectorDimArray}(undef,size(tlist)),Ti(tlist))
    
    C[1] = Ci + Cf
    
    # % Compute solution.
    for tt = 2:length(tlist)
        ti = tlist[tt-1]
        tf = tlist[tt]
        Ci = timestep_initial_condition(C[tt-1], μ, V, ti, tf)

        # Forcing contribution
        Cf = integrate_forcing( ti, tf, μ, V, B, source_history)

        # total
        C[tt] = Ci + Cf
    end # tt
    return real.(C) # Cut imaginary part which is zero to machine precision.
end

"""
    timestep_initial_condition(C, μ, V, ti, tf)

# Arguments
- `C::DimArray`: tracer distribution at `ti`
- `μ`: eigenvalue diagonal matrix
- `V`: eigenvector matrix
- `ti`: initial time
- `tf`: final time
# Returns
- `Cf::DimArray`: tracer distribution at `tf`
"""
timestep_initial_condition(C, μ, V, ti, tf) = real.( V * exp(Diagonal(μ)*(tf-ti)) / V * C )

"""
    forcing_integrand(t, tf, μ, V, B, source_history)

Integrand for boundary condition term in equation 10 (Haine et al., 2024).

# Arguments
- `t`: time
- `tf`: final time 
- `μ`: eigenvalue diagonal matrix
- `V`: eigenvector matrix
- `B`: boundary condition matrix
- `source_history::Function`: returns Dirichlet boundary condition at a given time
"""
forcing_integrand(t, tf, μ, V, B, source_history) = real.( V * exp(Diagonal(μ)*(tf-t)) / V * B * source_history(t))
    
"""
    integrate_forcing(t0, tf, μ, V, B, source_history)

Integrate boundary condition term in equation 10 (Haine et al., 2024).

# Arguments
- `t0`: initial time
- `tf`: final time 
- `μ`: eigenvalue diagonal matrix
- `V`: eigenvector matrix
- `B`: boundary condition matrix
- `source_history::Function`: returns Dirichlet boundary condition at a given time
"""
function integrate_forcing(t0, tf, μ, V, B, source_history)
    forcing_func(t) = forcing_integrand(t, tf, μ, V, B, source_history)

    # MATLAB: integral(integrand,ti,tf,'ArrayValued',true)
    integral, err = quadgk(forcing_func, t0, tf)
    (err < 1e-5) ? (return integral) : error("integration error too large")
end

"""
    tracer_timeseries(tracername, A, B, tlist, mbox1, vbox1; BD=nothing, halflife=nothing)

Simulate tracers and return tracer timeseries from one box.

# Arguments
- `tracername`: name of tracer
- `A`: tracer transport matrix
- `B`: boundary condition matrix
- `tlist`: list of times to save tracer concentration
- `mbox`: name of meridional box of interest
- `vbox`: name of vertical box of interest
- `BD=nothing`: Dirichlet boundary condition
- `halflife=nothing`: radioactive half life
"""
function tracer_timeseries(tracername, A, B, tlist, mbox1, vbox1; BD=nothing, halflife=nothing)

    if isnothing(halflife) && !isnothing(BD)
        return transient_tracer_timeseries(tracername, A, B, BD, tlist, mbox1, vbox1)
    elseif tracername == :argon39
        return steady_tracer_timeseries(tracername, A, B, halflife, tlist, mbox1, vbox1)
    elseif tracername == :iodine129 && !isnothing(BD)
        return transient_tracer_timeseries(tracername, A, B, BD, tlist, mbox1, vbox1, halflife = halflife)
    end
end

"""
    transient_tracer_timeseries(tracername, A, B, BD, tlist, mbox1, vbox1; halflife = nothing)

Simulate transient tracers and return tracer timeseries from one box.

# Arguments
- `tracername`: name of tracer
- `A`: tracer transport matrix
- `B`: boundary condition matrix
- `BD`: Dirichlet boundary condition
- `tlist`: list of times to save tracer concentration
- `mbox`: name of meridional box of interest
- `vbox`: name of vertical box of interest
- `halflife=nothing`: radioactive half life
"""
function transient_tracer_timeseries(tracername, A, B, BD, tlist, mbox1, vbox1; halflife = nothing)

    # fixed parameters for transient tracers
    if tracername == :iodine129
        box2_box1_ratio = 0.25
    elseif  (tracername == :CFC11NH) ||
        (tracername == :CFC12NH) ||
        (tracername == :SF6NH)
        box2_box1_ratio = 0.75
    else
        error("transient tracer not implemented")
    end

    # all tracers start with zero boundary conditions
    C₀ = zeros(model_dimensions())
	
    source_history_func(t) =  tracer_source_history(t,
	tracername,
	box2_box1_ratio,
	BD,
    )
    
    Cevolve = evolve_concentration(C₀, 
	A,
	B,
	tlist, 
	source_history_func,
        halflife = halflife)
    	
    return [Cevolve[t][At(mbox1),At(vbox1)] for t in eachindex(tlist)]

end

"""
    steady_tracer_timeseries(tracername, A, B, halflife, tlist, mbox1, vbox1)

Simulate non-transient tracers and return tracer timeseries from one box.

# Arguments
- `tracername`: name of tracer
- `A`: tracer transport matrix
- `B`: boundary condition matrix
- `halflife`: radioactive half life
- `BD`: Dirichlet boundary condition
- `tlist`: list of times to save tracer concentration
- `mbox`: name of meridional box of interest
- `vbox`: name of vertical box of interest
"""
function steady_tracer_timeseries(tracername, A, B, halflife, tlist, mbox1, vbox1)

    C₀ = ones(model_dimensions()) # initial conditions: faster spinup

    if tracername == :argon39
        box2_box1_ratio = 1 
    
        source_history_func(t) =  tracer_source_history(t,
	    tracername,
	    box2_box1_ratio
        )

        Cevolve = evolve_concentration(C₀, 
	    A,
	    B,
	    tlist, 
	    source_history_func;
	    halflife = halflife)

    else
        error("only implemented for argon-39")
    end
    
    return [Cevolve[t][At(mbox1),At(vbox1)] for t in eachindex(tlist)]

end


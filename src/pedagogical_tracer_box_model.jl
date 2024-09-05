# follows the design of Tom Haine's pedagaogical box model

kg = u"kg"
m  = u"m"
yr = u"yr"
Tg = u"Tg"
s  = u"s"
pmol = u"pmol"
fmol = u"fmol"
nmol = u"nmol"

@dim Tracer "tracer"
@dim Meridional "meridional location"
@dim Vertical "vertical location"
@dim Global "global quantity"

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
model_dimensions() = (Meridional(meridional_names()),Vertical(vertical_names())) 

function boundary_dimensions()
    meridional_boundary = meridional_names()[1:2]
    vertical_boundary = [vertical_names()[1]] # add brackets to keep as vector
    return  (Meridional(meridional_boundary), Vertical(vertical_boundary))
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

allequal(x) = all(y -> y == first(x), x)

maximum_timescale(μ) = -1/real(last(last(μ)))

function watermass_fraction(μ, V, B; alg=:forward)
    if alg == :forward
        return watermass_fraction_forward(μ, V, B)
    elseif alg == :adjoint 
        return watermass_fraction_adjoint(μ, V, B)
    elseif alg == :residence
        return watermass_fraction_residence(μ, V, B)
    else
        error("not yet implemented")
    end
end

watermass_fraction_forward(μ, V, B) = - real.(V / μ / V * B)

watermass_fraction_adjoint(μ, V, B) = - real.(transpose(B) * V / μ / V)

function watermass_fraction_residence(μ, V, B)
    # real(    B'*V/(D.^2)/V*B)
    μ_diag = diag(μ)
    μ2_diag = μ_diag.^2
    μ2 = DiagonalDimArray(μ2_diag,dims(μ))

    Nb = size(V) # number of boxes
    return real.( transpose(B) * V / μ2 / V * B) ./ Nb
end

function mean_age(μ, V, B; alg=:forward)
    if alg == :forward
        return mean_age_forward(μ, V, B)
    elseif alg == :adjoint 
        return mean_age_adjoint(μ, V, B)
    elseif alg == :residence
        return mean_age_residence(μ, V, B)
    else
        error("not yet implemented")
    end
end

function mean_age_forward(μ, V, B)
    μ_diag = diag(μ)
    μ2_diag = μ_diag.^2
    μ2 = DiagonalDimArray(μ2_diag,dims(μ))

    # use  real to get rid of very small complex parts
    # ideally, would check that complex parts are small
    boundary_dims = dims(B)
    return real.(V / μ2 / V * B) * ones(boundary_dims)
end

function mean_age_adjoint(μ, V, B)
    # MATLAB: [1, 1]*real(    B'*V/(D.^2)/V)
    μ_diag = diag(μ)
    μ2_diag = μ_diag.^2
    μ2 = DiagonalDimArray(μ2_diag,dims(μ))

    # use a 1 x 2 matrix to avoid ambiguity with transpose operator
    ones_row_vector = MultipliableDimArray(ones(1,2),Global(["mean age"]),dims(B))
    
    a_tmp = ones_row_vector * real.(transpose(B) * V / μ2 / V) 

    # undo the extra complication of a Global dimension
    return DimArray(reshape(transpose(Matrix(a_tmp)),size(a_tmp)),dims(a_tmp))
end
#function mean_age_adjoint(A, B)
    # previous working method 
    #μ, V = eigen(transpose(A))
    #return mean_age(μ, V, B)
#end

function mean_age_residence(μ, V, B)
    # MATLAB: [1, 1]*real(-2.*B'*V/(D.^3)/V*B)*[1; 1]./boxModel.no_boxes
    μ_diag = diag(μ)
    μ3_diag = μ_diag.^3 
    μ3 = DiagonalDimArray(μ3_diag,dims(μ))

    # use a 1 x 2 matrix to avoid ambiguity with transpose operator
    ones_row_vector = MultipliableDimArray(ones(1,2),Global(["mean age"]),dims(B))

    tmp = -2 .* ones_row_vector * real.(transpose(B) * V / μ3 / V * B) * transpose(ones_row_vector) 

    Nb = size(V) # number of boxes
    return tmp ./ Nb
end

function ttd_width(μ, V, B; alg=:forward)
    if alg == :forward
        return ttd_width_forward(μ, V, B)
    elseif alg == :adjoint 
        return ttd_width_adjoint(μ, V, B)
    elseif alg == :residence
        return ttd_width_residence(μ, V, B)
    else
        error("not yet implemented")
    end
end

function ttd_width_forward(μ, V, B)

    # MATLAB: sqrt((real(-2.*V/(D.^3)/V*B)*[1; 1] - (Solution.fwd_mean_ages).^2)./2) ;
    μ_diag = diag(μ)
    μ3_diag = μ_diag.^3 
    μ3 = DiagonalDimArray(μ3_diag,dims(μ))
    
    Δ² =  -real.(V / μ3 / V * B) * ones(dims(B))
    Γ = mean_age(μ, V, B, alg=:forward)
    Δ² -= ((1//2) .* Γ.^2)
    return .√(Δ²)
end

function ttd_width_adjoint(μ, V, B)
    # MATLAB: sqrt(([1, 1]*real(-2.*B'*V/(D.^3)/V) - (Solution.adj_mean_ages).^2)./2)
    μ_diag = diag(μ)
    μ3_diag = μ_diag.^3 
    μ3 = DiagonalDimArray(μ3_diag,dims(μ))

    # use a 1 x 2 matrix to avoid ambiguity with transpose operator
    ones_row_vector = MultipliableDimArray(ones(1,2),Global(["mean age"]),dims(B))
    
    Δ_tmp = -2 .* ones_row_vector * real.(transpose(B) * V / μ3 / V) 

    Δ2 =  DimArray(reshape(transpose(Matrix(Δ_tmp)),size(Δ_tmp)),dims(Δ_tmp))

    Γ = mean_age(μ, V, B, alg=:adjoint)
    Δ2 .-= Γ.^2 
   
    Δ² = (1//2) .* Δ2
    return .√(Δ²)
end

function ttd_width_residence(μ, V, B)
# MATLAB: sqrt(([1, 1]*real( 6.*B'*V/(D.^4)/V*B)*[1; 1]./boxModel.no_boxes - Solution.RTD_mean_rt^2)/2) ;

    μ_diag = diag(μ)
    μ4_diag = μ_diag.^4 
    μ4 = DiagonalDimArray(μ4_diag,dims(μ))

    # use a 1 x 2 matrix to avoid ambiguity with transpose operator
    ones_row_vector = MultipliableDimArray(ones(1,2),Global(["mean age"]),dims(B))

    Nb = size(V) # number of boxes
    tmp = (6 ./ Nb) .* ones_row_vector * real.(transpose(B) * V / μ4 / V * B) * transpose(ones_row_vector) 
    Γ = mean_age(μ, V, B, alg=:residence)

    # get rid of DimArrays for this global quantity (think about improving code design here)
    return .√((1//2) .* (first(first(tmp)) - first(first(Γ))^2 ))
end

normalized_exponential_decay(t,Tmax) = (1/Tmax)*exp(-(t/Tmax))

# older working method
# function adjoint_ttd_width(A, B)
#     μ, V = eigen(transpose(A))
#     return ttd_width(μ, V, B)
# end

location_tracer_histories() = "https://github.com/ThomasHaine/Pedagogical-Tracer-Box-Model/raw/main/MATLAB/tracer_gas_histories.mat"
#https://github.com/ThomasHaine/Pedagogical-Tracer-Box-Model/raw/main/MATLAB/tracer_gas_histories.mat
#download_tracer_histories() = Downloads.Download(url_tracer_histories())

function read_tracer_histories()

    # download tracer history input (make this lazy)
    url = OceanGreensFunctionMethods.location_tracer_histories()
    !isdir(datadir()) && mkpath(datadir())
    matfile = Downloads.download(url,datadir("tracer_histories.mat"))

    file = matopen(matfile)

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
    :SF6NH => NoUnits,
    :SF6SH => NoUnits,
    :N2ONH => nmol/kg,
    :N2OSH => nmol/kg
    )
    
function tracer_point_source_history(tracername, BD)

    tracer_timeseries = BD[Tracer=At(tracername)] * tracer_units()[tracername]

    return linear_interpolation(
        first(DimensionalData.index(dims(tracer_timeseries))),
        tracer_timeseries)
end

function tracer_source_history(t, tracername, BD, box2_box1_ratio)

    source_func = tracer_point_source_history(tracername, BD)
    box1 = source_func(t)
    box2 = box2_box1_ratio * box1
    
    # replace this section with a function call.
    boundary_dims = boundary_dimensions()
    return DimArray(hcat([box1,box2]),boundary_dims)
end

function evolve_concentration(C₀, A, B, tlist, source_history; halflife = nothing)
# % Integrate forcing vector over time to compute the concentration history.
# % Find propagator by analytical expression using eigen-methods.

    μ, V = eigen(A)

    # initial condition contribution
    Ci = deepcopy(C₀)

    # forcing contribution
    Cf = zeros(dims(C₀))

    # total
    C = DimArray(Array{DimArray}(undef,size(tlist)),Ti(tlist))
    
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
    return real.(C) #real.(Ci + Cf)  # Cut imaginary part which is zero to machine precision.
end

# MATLAB: concsI(:,tt) =      V*expm(D.*(tf-ti))/V*(concsI(:,tt-1) + concsF(:,tt-1)) ;    % Initial condition contribution
function timestep_initial_condition(C, μ, V, ti, tf)

    matexp = MultipliableDimArray( exp(Matrix(μ*(tf-ti))), dims(μ), dims(μ))
    return real.( V * (matexp * (V\C))) # matlab code has right divide (?)
end

# MATLAB: integrand    = @(t) V*expm(D.*(tf-t ))/V*B*source_history(t) ;
function forcing_integrand(t, tf, μ, V, B, source_history)

    matexp = MultipliableDimArray( exp(Matrix(μ*(tf-t))), dims(μ), dims(μ))

    # annoying finding: parentheses matter in next line
    return real.( V * (matexp * (V \ (B*source_history(t))))) 
end

function integrate_forcing(t0, tf, μ, V, B, source_history)

    Bunit = unit(first(first(B)))
    forcing_func(t) = Bunit * forcing_integrand(t, tf, μ, V, ustrip.(B), source_history)

    # MATLAB: integral(integrand,ti,tf,'ArrayValued',true)
    integral, err = quadgk(forcing_func, t0, tf)

    (err < 1e-5) ? (return integral) : error("integration error too large")
end

function transient_tracer_timeseries(tracername, BD, A, B, tlist, mbox1, vbox1)

    # fixed parameters for transient tracers
    box2_box1_ratio = 0.75
    C₀ = zeros(model_dimensions())
	
    source_history_func(t) =  tracer_source_history(t,
	tracername,
	BD,
	box2_box1_ratio)

    Cevolve = evolve_concentration(C₀, 
	A,
	B,
	tlist, 
	source_history_func;
	halflife = nothing)
	
    return [Cevolve[t][At(mbox1),At(vbox1)] for t in eachindex(tlist)]

end

greens_function(t,A::DimMatrix{DM}) where DM <: DimMatrix{Q} where Q <: Quantity = exp(A*t)

function boundary_propagator(t,A::DimMatrix{DM}, B::DimMatrix{DM}; alg=:forward) where DM <: DimMatrix
if alg == :forward 
    return boundary_propagator_forward(t, A, B)
elseif alg == :adjoint
    return boundary_propagator_adjoint(t, A, B)
end
    error("boundary propagator method not implemented")
end

boundary_propagator_forward(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = greens_function(t,A)*B
boundary_propagator_adjoint(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = transpose(B)*greens_function(t,A)

function global_ttd(t, A::DimMatrix{DM}, B::DimMatrix{DM}; alg=:forward) where DM <: DimMatrix
    if alg == :forward 
        return global_ttd_forward(t, A, B)
    elseif alg == :adjoint
        return global_ttd_adjoint(t, A, B)
    else
        error("global ttd method not implemented")
    end
end

global_ttd_forward(t, A::DimMatrix{DM}, B::DimMatrix{DM}) where DM <: DimMatrix = greens_function(t,A)*B*ones(dims(B))
global_ttd_adjoint(t, A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = transpose(ones(dims(B))) *  boundary_propagator_adjoint(t,A,B)

# not normalized by number of boxes: consistent with manuscript?
residence_time(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = t * transpose(B)*greens_function(t,A)*B

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
function global_ttd_adjoint(t, A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix
    ones_row_vector = MultipliableDimArray(ones(1,2),Global(["mean age"]),dims(B))
    tmp = ones_row_vector *  boundary_propagator_adjoint(t,A,B)

    # undo the extra complication of a Global dimension
    return DimArray(reshape(transpose(Matrix(tmp)),size(tmp)),dims(tmp))
end

# not normalized by number of boxes: consistent with manuscript?
residence_time(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = t * transpose(B)*greens_function(t,A)*B


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
    # MATLAB: real(    B'*V/(D.^2)/V*B)
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

    Nb = length(V) # number of boxes

    # get rid of DimArrays for this global quantity (think about improving code design here)
    return first(first(tmp)) / Nb
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

    return .√((1//2) .* (first(first(tmp)) - Γ^2 ))
end

normalized_exponential_decay(t,Tmax) = (1/Tmax)*exp(-(t/Tmax))

# path density for a given time and location 
function path_density(μ, V, B, t, mbox, vbox)
    Φ(τ) = phi_function(μ, τ) # a useful closure
    D_mat = MultipliableDimArray(zeros(length(V), length(V)),model_dimensions(),model_dimensions())
    D_mat[At(mbox),At(vbox)][At(mbox),At(vbox)] = 1 
    D_mat_overline = V \ D_mat * V

    # need to define element-by-element product for MultipliableDimArrays
    elemental_product = MultipliableDimArray(Matrix(D_mat_overline).*Matrix(Φ(t)),
        dims(D_mat_overline), dims(D_mat_overline))

    #return real.( transpose(B) * V * (D_mat_overline .* Φ(t)) / V * B)
    return real.( transpose(B) * V * elemental_product / V * B)
end

function phi_function(μ, t)
    N = length(μ)
    eigen_dims = MultipliableDimArrays.Eigenmode(1:N)
    ϕ = MultipliableDimArray(zeros(ComplexF64, N, N)yr, eigen_dims, eigen_dims)
    μvals = diag(μ)
    for rr in 1:N
        for cc in 1:N
            μ_rr = μvals[rr]
            μ_cc = μvals[cc]
            if μ_rr ≠ μ_cc
                ϕ[cc][rr] = (exp(μ_cc*t) - exp(μ_rr*t))/(μ_cc - μ_rr)
            else
                ϕ[cc][rr] = t*exp(μ_rr*t)
            end
        end # cc
    end # rr
    return ϕ
end

function ideal_age(A, B; alg= :forward)
    if alg == :forward
        return ideal_age_forward(A, B)
    elseif alg == :adjoint 
        return ideal_age_adjoint(A, B)
    else
        error("not yet implemented")
    end
end

ideal_age_forward(A, B) = - A \ (B*zeros(boundary_dimensions())yr + ones(model_dimensions()))

# doesn't work due to conflict with DimensionalData.transpose that I don't want to overload
ideal_age_adjoint(A, B) = - transpose(A) \ (B*zeros(boundary_dimensions())yr + ones(model_dimensions()))

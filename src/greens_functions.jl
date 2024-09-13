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
    #ϕ = zeros(ComplexF64,model_dimensions())yr # \phi + TAB
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

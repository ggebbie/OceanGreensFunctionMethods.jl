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

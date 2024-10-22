"""
    greens_function(τ,A)

Green's function for a box model (for steady transport given by the matrix 𝐀 for response at time t to a source at time t′ where τ = t - t′): the matrix exponential function of the elapsed time between the source time and field time:
```math
{\\bf G}(\\tau) = e^{ {\\bf A} \\tau}
```
where 𝐆(t) is a  N × N matrix with the spatial locations of field points (boxes) down its N rows and source points (boxes) along its N columns. Thus, the element 𝐆{i,j}(τ) quantifies transfer from a source at time t′ in box j to receiver at time t in box i.
"""
greens_function(τ,A::AbstractMatrix) = exp(A*τ)
# is specialized code (below) necessary?
#greens_function(τ,A::DimMatrix{DM}) where DM <: DimMatrix{Q} where Q <: Quantity = exp(A*τ)

"""
    boundary_propagator(τ, A, B; alg=:forward)

Forward and adjoint boundary propagators.

# Forward boundary propagator

    boundary_propagator(τ, A, B, alg=:forward)

The (forward) boundary propagator is the box model surface-to-interior transit time distribution (TTD) over transit times τ = t - t′, as given by equation 88 of Haine et al. (2024):
```math
{\\bf G}' (\\tau) = {\\bf G} (\\tau) ~ {\\bf B}
```
The N × Nₛ 𝐆′(τ) matrix quantifies transfer from the Nₛ components of the surface forcing to the N boxes with transit time τ.

# Adjoint boundary propagator

    boundary_propagator(τ, A, B, alg=:adjoint)

The box model adjoint boundary propagator (interior-to-surface TTD over transit time τ† = t″ - t, where t″ ≥ t is the time of the adjoint source; equation 93 of Haine et al., 2024) is
```math
{\\bf G}'^{\\dagger} (\\tau^{\\dagger} )  = {\\bf B}^{T}~ {\\bf G} (\\tau^{\\dagger}).
```
This Nₛ × N 𝐆′†(τ†) matrix quantifies transfer from the N interior boxes to the Nₛ surface boxes with transit time τ†.
"""
function boundary_propagator(τ, A::AbstractMatrix, B::AbstractMatrix; alg=:forward) 
#function boundary_propagator(τ, A::DimMatrix{DM}, B::DimMatrix{DM}; alg=:forward) where DM <: DimMatrix
if alg == :forward 
    return boundary_propagator_forward(τ, A, B)
elseif alg == :adjoint
    return boundary_propagator_adjoint(τ, A, B)
end
    error("boundary propagator method not implemented")
end

"""
    boundary_propagator_forward(t,A,B)
"""
boundary_propagator_forward(t,A::AbstractMatrix, B::AbstractMatrix) = greens_function(t,A)*B
#boundary_propagator_forward(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = greens_function(t,A)*B

"""
    boundary_propagator_adjoint(t,A,B)
"""
boundary_propagator_adjoint(t, A::AbstractMatrix, B::AbstractMatrix) = transpose(B)*greens_function(t,A)
#boundary_propagator_adjoint(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = transpose(B)*greens_function(t,A)

"""
    global_ttd(t, A, B; alg=:forward)

Forward and adjoint global transit time distributions (TTDs).

# Forward Global TTD

The forward global (total) TTD is the sum of surface-to-interior TTDs (equation 90 of Haine et al., 2024): 
```math
{\\cal G} (t) = {\\bf G} (t) ~ {\\bf B} ~ {\\bf 1}_{N_S},
```
where the product with the Ns × 1 column vector of ones (i.e., last matrix in previous equation) computes the sum over surface boxes. This expression yields an N × 1 column vector that is normalized for each box.

# Adjoint Global TTD

The adjoint global (total) TTD is the sum of interior-to-surface TTDs. 
"""
function global_ttd(t, A::AbstractMatrix, B::AbstractMatrix; alg=:forward) 
#function global_ttd(t, A::DimMatrix{DM}, B::DimMatrix{DM}; alg=:forward) where DM <: DimMatrix
    if alg == :forward 
        return global_ttd_forward(t, A, B)
    elseif alg == :adjoint
        return global_ttd_adjoint(t, A, B)
    else
        error("global ttd method not implemented")
    end
end

"""
    global_ttd_forward(t, A, B)
"""
global_ttd_forward(t, A::AbstractMatrix, B::AbstractMatrix) = greens_function(t,A)*B*ones(dims(B))
#global_ttd_forward(t, A::DimMatrix{DM}, B::DimMatrix{DM}) where DM <: DimMatrix = greens_function(t,A)*B*ones(dims(B))

"""
    global_ttd_adjoint(t, A, B)
"""
function global_ttd_adjoint(t, A::AbstractMatrix, B::AbstractMatrix)
    #function global_ttd_adjoint(t, A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix
    boundary_dims = domainsize(B)

    # ones_row_vector = AlgebraicArray(ones(1,2),Global(["mean age"]),dims(B))
    # tmp = ones_row_vector *  boundary_propagator_adjoint(t,A,B)
    return transpose( transpose(ones(boundary_dims, :VectorArray)) * boundary_propagator_adjoint(t,A,B) )
    
    # undo the extra complication of a Global dimension
    #return AlgebraicArray(transpose(Matrix(tmp)),dims(tmp))
    #return VectorArray(DimArray(reshape(transpose(Matrix(tmp)),size(tmp)),dims(tmp)))
    #return transpose(tmp)
end

"""
    residence_time(t, A, B)

The surface-to-surface residence-time distribution (RTD) is (equations 94 and 95 of Haine et al., 2024):
```math
{\\bf R} (\\tau) = 
\\frac{1}{N} \\int_{t - \\tau}^{t} {\\bf G}'^{\\dagger} (t^* + \\tau - t)  ~ {\\bf G}' (t - t^*) ~ d t ^*
```
or
```math
{\\bf R} (\\tau)  = \\frac{\\tau}{N} {\\bf B}^{T}  ~ {\\bf G} (\\tau) ~ {\\bf B},
```
where N is the number of boxes, G(τ) is the forward Green's function and B is the boundary matrix. 
The Ns × Ns R(τ) matrix quantifies transfer from the Ns surface boxes back to the Ns surface boxes with residence time τ (element R{i,j}(τ) quantifies transfer from entry box j to exit box i).

Note: not normalized by number of boxes in this code: consistent with manuscript?
"""
residence_time(t, A::AbstractMatrix, B::AbstractMatrix) = t*transpose(B)*greens_function(t,A)*B
#residence_time(t,A::DimMatrix{DM},B::DimMatrix{DM}) where DM <: DimMatrix = t * transpose(B)*greens_function(t,A)*B

"""
    maximum_timescale(μ)

Return `Tmax` for the eigenvalues μ. The matrix exponential of 𝐀τ has asymptotic properties because G(t) must eventually decay exponentially with timescale 
```math
T_{max} = -1/\\mu_{min},        
```
where μmin is the eigenvalue with smallest real part. Thus, the Green's function has a maximum timescale of Tmax which is larger than all other transport timescales.
"""
maximum_timescale(μ) = -1/real(last(μ))

"""
    watermass_fraction(μ, V, B; alg=:forward)

Forward, adjoint, and residence-time water-mass fractions.

# Forward water-mass fraction

    watermass_fraction(μ, V, B, alg=:forward)

The water mass fractions are (equation 89 of Haine et al., 2024)
```math
{\\bf a}  = \\int_0^{\\infty} {\\bf G} (\\tau) ~ {\\bf B} ~ d \\tau
```
or
```math
{\\bf a} = -{\\bf V} ~ \\mu^{-1} ~ {\\bf V}^{-1} ~ {\\bf B} , 
```
which is an N × Ns matrix with the interior boxes down the rows and the surface sources across the rows.

# Adjoint water-mass fraction

    watermass_fraction(μ, V, B, alg=:adjoint)

Fraction of water that will return to the surface in a particular box.

# Residence-time water-mass fraction

    watermass_fraction(μ, V, B, alg=:residence)

Fraction of water that leaves a particular box and returns in another box.
"""
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

"""
    watermass_fraction_forward(μ, V, B)
"""
watermass_fraction_forward(μ, V, B) = - real(V/ Diagonal(μ) / V * B)

"""
    watermass_fraction_adjoint(μ, V, B)
"""
watermass_fraction_adjoint(μ, V, B) = - real(transpose(B) * V / Diagonal(μ) / V)

"""
    watermass_fraction_residence(μ, V, B)
"""
function watermass_fraction_residence(μ, V, B)
    # MATLAB: real(    B'*V/(D.^2)/V*B)
    D2 = Diagonal(μ.^2)
    Nb = length(V) # number of boxes
    return real( transpose(B) * V / D2 / V * B) / Nb
end

"""
    mean_age(μ, V, B; alg=:forward)

Mean age of the forward TTDs, adjoint TTDs, and residence-time distributions.
# Arguments
- `μ`: eigenvalues vector
- `V`: eigenvector matrix
- `B`: boundary matrix
- `alg=:forward`: algorithm (optional)

# Forward mean age

    mean_age(μ, V, B, alg=:forward)

The mean transit time 𝚪 (mean age) is (equation 92 of Haine et al., 2004),
```math
{\\bf \\Gamma} = {\\bf V} ~ \\mu^{-2} ~ {\\bf V}^{-1} ~ {\\bf B} ~ {\\bf 1}_{N_S},
```
which is an N × 1 vector for each box (and which also equals the ideal age).

# Adjoint mean age

    mean_age(μ, V, B, alg=:adjoint)

# Residence-time mean age

    mean_age(μ, V, B, alg=:residence)
"""
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

"""
    mean_age_forward(μ, V, B)
"""
function mean_age_forward(μ, V, B)
    # μ_diag = diag(μ)
    # μ2_diag = μ_diag.^2
    # μ2 = DiagonalDimArray(μ2_diag,dims(μ))
    D2 = Diagonal(μ.^2)
    
    # use  real to get rid of very small complex parts
    # ideally, would check that complex parts are small
    boundary_dims = domainsize(B)
    return real(V / D2 / V ) * B * ones(boundary_dims, :VectorArray)
end

"""
    mean_age_adjoint(μ, V, B)
"""
function mean_age_adjoint(μ, V, B)
    # MATLAB: [1, 1]*real(    B'*V/(D.^2)/V)
    D2 = Diagonal(μ.^2)

    # use a 1 x 2 matrix to avoid ambiguity with transpose operator
    # ones_row_vector = AlgebraicArray(ones(1,2),Global(["mean age"]),dims(B))
    # a_tmp = ones_row_vector * real(transpose(B) * V / D2 / V)
    boundary_dims = domainsize(B)
    Γ = transpose(ones(boundary_dims, :VectorArray))  * real(transpose(B) * V / D2 / V) 

    # undo the extra complication of a Global dimension
    return transpose(Γ)
end

"""
    mean_age_residence(μ, V, B)
"""
function mean_age_residence(μ, V, B)
    # MATLAB: [1, 1]*real(-2.*B'*V/(D.^3)/V*B)*[1; 1]./boxModel.no_boxes
    D3 = Diagonal(μ.^3)
    boundary_dims = domainsize(B)

    Γ = -2 * transpose(ones(boundary_dims, :VectorArray))*
        real(transpose(B) * V / D3 / V * B) *
        ones(boundary_dims, :VectorArray)

    Nb = length(V) # number of boxes
    return Γ / Nb
end

"""
    ttd_width(μ, V, B; alg=:forward)

Width of the forward TTDs, adjoint TTDs, and residence-time distributions.

# Arguments
- `μ`: eigenvalue diagonal matrix
- `V`: eigenvector matrix
- `B`: boundary matrix
- `alg=:forward`: algorithm (optional)
# Returns
- `Δ`: TTD width

# Width of forward TTD

    ttd_width(μ, V, B, alg=:forward)

The TTD width is given by (equation 92 of Haine et al., 2024),
```math
2 {\\bf \\Delta}^2  = -2 ~ {\\bf V} ~ \\mu^{-3} ~ {\\bf V}^{-1} ~ {\\bf B} ~ {\\bf 1}_{N_S}  - {\\bf \\Gamma}^2,
```
which is a N × 1 vector for each box.

# Adjoint mean age

    mean_age(μ, V, B, alg=:adjoint)

# Residence-time mean age

    mean_age(μ, V, B, alg=:residence)
"""
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

"""
    ttd_width_forward(μ, V, B)
"""
function ttd_width_forward(μ, V, B)
    # MATLAB: sqrt((real(-2.*V/(D.^3)/V*B)*[1; 1] - (Solution.fwd_mean_ages).^2)./2) ;
    D3 = Diagonal(μ.^3)

    Δ² =  - real(V / D3 / V * B) * ones(domainsize(B), :VectorArray)
    Γ = mean_age(μ, V, B, alg=:forward)
    Δ² -= ((1//2) .* Γ.^2)
    return .√(Δ²)
end

"""
    ttd_width_adjoint(μ, V, B)
"""
function ttd_width_adjoint(μ, V, B)
    # MATLAB: sqrt(([1, 1]*real(-2.*B'*V/(D.^3)/V) - (Solution.adj_mean_ages).^2)./2)
    D3 = Diagonal(μ.^3)
    boundary_dims = domainsize(B)
    Δ = -2 * transpose( transpose(ones(boundary_dims, :VectorArray)) * real(transpose(B) * V / D3 / V) )
    Γ = mean_age(μ, V, B, alg=:adjoint)
    Δ .-= Γ.^2 
    Δ .*= (1//2) 
    return .√(Δ)
end

"""
    ttd_width_residence(μ, V, B)
"""
function ttd_width_residence(μ, V, B)
    # MATLAB: sqrt(([1, 1]*real( 6.*B'*V/(D.^4)/V*B)*[1; 1]./boxModel.no_boxes - Solution.RTD_mean_rt^2)/2) ;
    D4 = Diagonal(μ.^4)
    boundary_dims = domainsize(B)
    Nb = length(V) # number of boxes
    Δ2 = (6 / Nb) *
        transpose(ones(boundary_dims, :VectorArray)) *
        real(transpose(B) * V / D4 / V * B) *
        ones(boundary_dims, :VectorArray)

    Γ = mean_age(μ, V, B, alg=:residence)
    return .√((1//2) .* (Δ2 - Γ^2 ))
end

"""
    normalized_exponential_decay(t,Tmax)
"""
normalized_exponential_decay(t,Tmax) = (1/Tmax)*exp(-(t/Tmax))

"""
    path_density(μ, V, B, t, mbox, vbox)

# Arguments
- `μ`: eigenvalue diagonal matrix
- `V`: eigenvector matrix
- `B`: boundary matrix
- `t`: time
- `mbox`: name of meridional box of interest
- `vbox`: name of vertical box of interest
# Returns
- `E`: path density

The path density 𝐄_i(τ) for i ∈ 1 ... N is (equation 96 of Haine et al., 2024):
```math
{\\bf E}_i (\\tau)  = 
\\frac{1}{N} \\int_{t - \\tau}^{t} {\\bf G}'^{\\dagger} (t^* + \\tau - t) ~ {\\bf D}_i  ~ {\\bf G}' (t - t^*) ~ d t ^* , 
```
where 𝐃i is the N × N matrix unit of zeros with a single one at the i-th row and i-th column.
Therefore, 
```math
{\\bf E}_i (\\tau)  = \\frac{1}{N} \\int_{0}^{\\tau} {\\bf G}'^{\\dagger} (t') ~ {\\bf D}_i  ~ {\\bf G}' (\\tau - t') ~ d t '
```
and
```math
{\\bf E}_i (\\tau) = \\frac{1}{N} {\\bf B}^{T} \\int_{0}^{\\tau} ~ e^{{\\bf A} t'} ~ {\\bf D}_i  e^{{\\bf A} (\tau - t')} ~ d t ' {\\bf B}
```
and
```math
{\\bf E}_i (\\tau) = \\frac{1}{N}{\\bf B}^{T} ~ {\\bf V} \\left( \\overline{\\bf D}_i \\circ \\Phi (t) \\right) {\\bf V}^{-1} ~ {\\bf B}
```
where ϕ is defined in equation 100 of Haine et al. (2024). For a particular interior box i, 𝐄_i(τ) is the density of pathways between all combinations of surface entry and surface exit boxes over total residence time τ.
"""
function path_density(μ, V, B, t, mbox, vbox)
    Φ(τ) = phi_function(τ, μ) # a useful closure
    D_mat = AlgebraicArray(zeros(length(V), length(V)),model_dimensions(),model_dimensions())
    D_mat[At(mbox),At(vbox)][At(mbox),At(vbox)] = 1 
    D_mat_overline = V \ D_mat * V

    # check for element-by-element product to simplify 
    elemental_product = AlgebraicArray(Matrix(D_mat_overline).*Matrix(Φ(t)),
        dims(D_mat_overline), dims(D_mat_overline))

    #return real.( transpose(B) * V * (D_mat_overline .* Φ(t)) / V * B)
    return real( transpose(B) * V * elemental_product / V * B)
end

"""
    phi_function(t, μ)
"""
function phi_function(t, μ)
    N = length(μ) # correct translation for eigenvalue vector?
    #N = (length(μ))^2 # correct translation for eigenvalue vector?
    #eigen_dims = AlgebraicArrays.Eigenmode(1:N)
    eigen_dims = Eigenmode(1:N)
    ϕ = AlgebraicArray(zeros(ComplexF64, N, N)yr, eigen_dims, eigen_dims)
    
    #μvals = diag(μ)
    for rr in 1:N
        println("rr",rr)
        for cc in 1:N
            println("rr",rr)
            μ_rr = μ[rr]
            μ_cc = μ[cc]
            # μ_rr = μvals[rr]
            # μ_cc = μvals[cc]
            if μ_rr ≠ μ_cc
                ϕ[cc][rr] = (exp(μ_cc*t) - exp(μ_rr*t))/(μ_cc - μ_rr)
            else
                ϕ[cc][rr] = t*exp(μ_rr*t)
            end
        end # cc
    end # rr
    return ϕ
end

"""
    ideal_age(A, B; alg= :forward)
"""
function ideal_age(A, B; alg= :forward)
    if alg == :forward
        return ideal_age_forward(A, B)
    elseif alg == :adjoint 
        return ideal_age_adjoint(A, B)
    else
        error("not yet implemented")
    end
end

"""
    ideal_age_forward(A, B)
"""
ideal_age_forward(A, B) = - A \ (B*zeros(domainsize(B),:VectorArray)yr + ones(domainsize(A),:VectorArray))

"""
    ideal_age_adjoint(A, B)
"""
ideal_age_adjoint(A, B) = - transpose(A) \ (B*zeros(domainsize(B), :VectorArray)yr + ones(domainsize(A), :VectorArray))

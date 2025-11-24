"""
Evaluate a polynomial at Θ, given coefficients for Θ¹, Θ², ..., Θ⁶.
Assumes coeffs = [a₁, a₂, ..., a₆] for Θ¹, Θ², ..., Θ⁶.
"""
function eval_poly_theta(Θ, coeffs)
    # coeffs[1]*Θ + coeffs[2]*Θ^2 + ... + coeffs[6]*Θ^6
    return sum(coeffs[i] * Θ^i for i in 1:length(coeffs))
end

"""
Evaluate the derivative of the polynomial at Θ.
Returns d/dΘ [coeffs[1]*Θ + coeffs[2]*Θ^2 + ... + coeffs[6]*Θ^6]
"""
function eval_poly_theta_deriv(Θ, coeffs)
    # coeffs[1]*1*Θ^0 + coeffs[2]*2*Θ^1 + ... + coeffs[6]*6*Θ^5
    return sum(coeffs[i] * i * Θ^(i-1) for i in 1:length(coeffs))
end

"""
Generic explicit RK interpolation.
Arguments:
- Θ: interpolation parameter (0 ≤ Θ ≤ 1)
- dt: time step
- y₀: initial value
- k: stage derivatives (array of vectors, one per stage)
- tableau: matrix of interpolation coefficients (nstages × ncoeffs)
  - Row i contains polynomial coefficients for stage i
  - Columns are coefficients for Θ⁰, Θ¹, Θ², ..., Θⁿ
- idxs: indices (optional, for partial interpolation)
- order: 0 for value, 1 for derivative
- k_indices: which stages from k to use (defaults to 1:nstages)
"""
function generic_interpolant(Θ, dt, y₀, k, tableau; idxs=nothing, order=0, k_indices=nothing)
    # Determine number of stages used in interpolation from tableau size
    nstages = size(tableau, 1)

    # Default: use all stages sequentially k[1], k[2], ..., k[nstages]
    if isnothing(k_indices)
        k_indices = 1:nstages
    end

    # For each stage, evaluate the polynomial or its derivative
    b = if order == 0
        [eval_poly_theta(Θ, tableau[i, :]) for i in 1:nstages]
    else
        [eval_poly_theta_deriv(Θ, tableau[i, :]) for i in 1:nstages]
    end

    # Compute the interpolation sum
    if isnothing(idxs)
        # Full vector
        interp_sum = sum(k[k_indices[i]] * b[i] for i in 1:nstages)
        if order == 0
            return y₀ + dt * interp_sum
        else
            return interp_sum
        end
    else
        # Indexed
        interp_sum = sum(k[k_indices[i]][idxs] * b[i] for i in 1:nstages)
        if order == 0
            return y₀[idxs] + dt * interp_sum
        else
            return interp_sum
        end
    end
end
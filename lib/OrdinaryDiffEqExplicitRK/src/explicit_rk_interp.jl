"""
Generic interpolation for Runge-Kutta methods.
Arguments:
- Θ: interpolation parameter (0 ≤ Θ ≤ 1)
- dt: time step
- y₀: initial value
- k: stage derivatives (vector of vectors, one per component)
- tableau: coefficient matrix where each row contains polynomial coefficients for a stage
          Each row i contains [a₀, a₁, a₂, ...] for polynomial aᵢ₀ + aᵢ₁*Θ + aᵢ₂*Θ² + ...
- idxs: indices (optional, for partial interpolation)
- order: 0 for value, 1 for derivative
"""
function generic_interpolant(Θ, dt, y₀, k, tableau; idxs=nothing, order=0)
    # Determine the number of stages based on the tableau size
    num_stages = size(tableau, 1)
    num_coeffs = size(tableau, 2)

    # For each stage, evaluate the polynomial or its derivative
    b = if order == 0
        # Use builtin evalpoly for polynomial evaluation: a₀ + a₁*Θ + a₂*Θ² + ...
        [@evalpoly(Θ, tableau[i,:]...) for i in 1:num_stages]
    else
        # For derivative: d/dΘ [a₀ + a₁*Θ + a₂*Θ² + ...] = a₁ + 2*a₂*Θ + 3*a₃*Θ² + ...
        [@evalpoly(Θ, [j * tableau[i, j+1] for j in 1:(num_coeffs-1)]...) for i in 1:num_stages]
    end

    # Compute the interpolation sum
    if isnothing(idxs)
        # Full vector
        interp_sum = sum(k[i] * b[i] for i in 1:num_stages)
        if order == 0
            return y₀ + dt * interp_sum
        else
            return interp_sum
        end
    else
        # Indexed
        interp_sum = sum(k[i][idxs] * b[i] for i in 1:num_stages)
        if order == 0
            return y₀[idxs] + dt * interp_sum
        else
            return interp_sum
        end
    end
end
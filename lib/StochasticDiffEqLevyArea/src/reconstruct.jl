"""
    reconstruct_path(dW, h, coeffs, t_points)

Reconstruct the Brownian motion at specified time points within [0, h] using the
Karhunen-Loève (Fourier) expansion:

    W_j(t) = (t/h)*dW_j + Σ_{k=1}^n a_{j,k} * √(2h)/(kπ) * sin(kπt/h)

where a_{j,k} are the Fourier coefficients stored in `coeffs`.

Returns a vector of m-dimensional vectors, one per time point.

Note: The coefficient layout differs by algorithm. For Fourier/Milstein/Wiktorsson,
`coeffs.X[k,j]` = a_{j,k}. For MronRoe, `coeffs.X[j,k]` = a_{j,k}.
"""
function reconstruct_path(
        dW::AbstractVector{T}, h::Real, coeffs::LevyAreaCoefficients{T},
        t_points::AbstractVector
    ) where {T <: AbstractFloat}
    m = length(dW)
    n = coeffs.n
    X = coeffs.X
    is_transposed = size(X, 1) == m && size(X, 2) == n  # MronRoe layout

    result = Vector{Vector{T}}(undef, length(t_points))
    for (idx, t) in enumerate(t_points)
        W_t = (t / h) .* dW
        for k in 1:n
            coeff = √(T(2) * h) / (k * T(π))
            s = sin(k * T(π) * t / h)
            for j in 1:m
                a_jk = is_transposed ? X[j, k] : X[k, j]
                W_t[j] += a_jk * coeff * s
            end
        end
        result[idx] = W_t
    end
    return result
end

"""
    iterated_integrals_subinterval(dW, h, coeffs, t_start, t_end;
        n_quadrature=64, ito_correction=true)

Compute iterated Stratonovich integrals over the sub-interval [t_start, t_end] ⊂ [0, h]
by reconstructing the Brownian path from Fourier coefficients and computing via
Riemann sums.

Returns an m×m matrix J where:
- With ito_correction=true: J_{jk} = Ito iterated integral ∫∫ dW_j dW_k
- With ito_correction=false: J_{jk} = Stratonovich iterated integral

The accuracy improves with `n_quadrature` (number of evaluation points within
the sub-interval) and with the truncation level `n` in the coefficients.
"""
function iterated_integrals_subinterval(
        dW::AbstractVector{T}, h::Real,
        coeffs::LevyAreaCoefficients{T},
        t_start::Real, t_end::Real;
        n_quadrature::Int = 64,
        ito_correction::Bool = true
    ) where {T <: AbstractFloat}
    m = length(dW)
    dt_sub = t_end - t_start

    # Generate quadrature points
    t_points = range(T(t_start), T(t_end), length = n_quadrature + 1)

    # Reconstruct path at quadrature points
    W_values = reconstruct_path(dW, h, coeffs, collect(t_points))

    # The reconstructed path is a truncated Fourier series (smooth function),
    # so the Riemann sum naturally gives the Stratonovich integral
    # (the quadratic variation of a smooth function is zero).
    # I_{jk}^{Strat} ≈ Σ_n dW_k^{(n)} * (W_j^{cumsum} + dW_j^{(n)}/2)
    # which equals Σ_n W_cumsum_j * dW_k^n + (1/2) * Σ_n dW_j^n * dW_k^n
    # The second term → (1/2)*δ_{jk}*dt for Brownian motion, but → 0 for smooth paths.
    # So the left-point Riemann sum gives the Itô integral + residual ≈ Stratonovich.

    # Compute via left-point Riemann sum (gives Stratonovich-like result for smooth paths):
    J = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for s in 2:length(W_values)
        dWs = W_values[s] .- W_values[s - 1]
        for k in 1:m
            for j in 1:m
                J[j, k] += W_cumsum[j] * dWs[k]
            end
        end
        W_cumsum .+= dWs
    end

    # The Riemann sum gives approximately: J_{jk} = I_{jk}^{Ito} + (1/2)*QV_{jk}
    # where QV_{jk} is the cross quadratic variation. For the truncated path,
    # QV → 0 (smooth), so J ≈ I_{Ito} with the truncation's Ito correction missing.
    # The true Ito integral has I_{jj} = (dW_j² - dt)/2 on diagonal.
    # Our Riemann sum gives J_{jj} ≈ dW_j²/2 (missing -dt/2).
    # To get the correct Stratonovich integral: J_{jk}^{Strat} = J_{jk} + (1/2)*δ_{jk}*dt_sub
    # To get the correct Ito integral: I_{jk} = J_{jk} - (1/2)*δ_{jk}*dt_sub
    # But our J already approximates the Ito integral (left-point sum without QV),
    # so to get Stratonovich we add dt/2, and J itself is approximately Ito.

    # Actually: for the truncated (smooth) path, the left-point Riemann sum equals
    # the midpoint sum equals the right-point sum (they all converge to the same
    # Riemann integral as N→∞). This Riemann integral is:
    # ∫_0^T W_j(t) dW_k(t) in the classical (non-stochastic) sense.
    # For the true BM: left-point → Itô, midpoint → Stratonovich.
    # For smooth approximation: they all → same value, which equals the Stratonovich
    # integral of the smooth path (= Itô integral + 0, since QV=0).

    # The off-diagonal elements are correct (no QV correction needed).
    # The diagonal needs: I_{jj}^{Ito} = ∫ W_j dW_j = (W_j(T)² - W_j(0)²)/2 - [W_j]/2
    # For the true BM: [W_j] = dt, so I_{jj} = (dW_j² - dt)/2
    # For smooth path: [W_j] = 0, so the integral = dW_j²/2

    # Correction: subtract dt_sub/2 from diagonal to match the true Itô integral
    if ito_correction
        for j in 1:m
            J[j, j] -= dt_sub / 2
        end
    end

    return J
end

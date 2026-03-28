"""
    _unpack_dZ_to_coefficients(dZ, m, alg_p)

Unpack the flat dZ vector (N(0,1) white noise) into LevyAreaCoefficients for MronRoe.
The dZ layout matches generate_coefficients(m, n, MronRoe(), rng): X (m×n), Y (n×m),
Ψ (m), tail ((m²-m)/2), packed contiguously.
"""
function _unpack_dZ_to_coefficients(dZ, m, alg_p)
    T = eltype(dZ)
    n_total = length(dZ)
    # MronRoe: norv = 2*m*n + (m²+m)/2
    # Solve for n: n = (n_total - (m²+m)/2) / (2m)
    n_tail = (m^2 + m) ÷ 2
    n = (n_total - n_tail) ÷ (2 * m)
    n > 0 || return nothing

    # Unpack: X is m×n, Y is n×m, tail is (m²+m)/2
    offset = 0
    X = reshape(dZ[(offset + 1):(offset + m * n)], m, n)
    offset += m * n
    Y = reshape(dZ[(offset + 1):(offset + n * m)], n, m)
    offset += n * m
    tail = dZ[(offset + 1):(offset + n_tail)]

    return LevyAreaCoefficients{T}(collect(X), collect(Y), collect(tail), m, n)
end

"""
    _compute_iterated_I(dt, dW, dZ, W_noise, cache, alg)

Compute Stratonovich iterated integrals J. Strategy:

1. **RSwM3 S₂ stack**: After step rejection, RSwM3 decomposes the current step's
   noise into sub-intervals stored in S₂. Compute iterated integrals from these
   exact sub-interval dW values. Handles multiple rejections naturally since
   `reject_step!` keeps S₂ updated.

2. **NoiseWrapper/NoiseGrid sub-grid**: If the noise process has fine-grid W values
   (from a prior solve), compute iterated integrals from sub-grid dW increments.
   This ensures consistency between fine and coarse solutions for convergence testing.

3. **Fresh step with dZ coefficients**: Unpack dZ into MronRoe Fourier coefficients,
   compute levyarea deterministically.

4. **Fallback**: Return nothing → legacy path for diagonal/scalar noise.
"""
function _compute_iterated_I(dt, dW, dZ, W_noise, alg)
    # Only for non-diagonal vector noise (dW must be a Vector, not Matrix or Number)
    dW isa AbstractVector || return nothing

    m = length(dW)

    # Strategy 1: RSwM3 S₂ stack (adaptive, post-rejection)
    J_s2 = _compute_II_from_S2(W_noise, m, dt)
    if J_s2 !== nothing
        return J_s2
    end

    # Strategy 2: NoiseWrapper/NoiseGrid sub-grid (convergence testing)
    J_subgrid = _compute_II_from_grid(W_noise, m, dt)
    if J_subgrid !== nothing
        return J_subgrid
    end

    # Check if we have usable dZ for coefficient-based computation
    if dZ === nothing || length(dZ) < 2 * m
        return nothing
    end

    # Strategy 3: Fresh step — unpack dZ into Fourier coefficients
    coeffs = _unpack_dZ_to_coefficients(dZ, m, alg.p)
    coeffs === nothing && return nothing

    # Compute full-step Lévy area from coefficients (Stratonovich)
    Wn = dW / √dt
    A = levyarea(Wn, coeffs.n, MronRoe(), coeffs)
    J = 1 // 2 * dW .* dW' .+ dt .* A
    return J
end

"""
Compute Stratonovich iterated integrals from RSwM3's S₂ stack.
S₂ holds the sub-interval decomposition (dt_i, dW_i, dZ_i) of the current
step's noise. After rejection, S₂ has ≥2 entries. With only 1 entry (fresh
step, no rejections), fall back to nothing so the dZ coefficient path is used.

The iterated integral over the composite step [0,h] = ∪_i [t_i, t_{i+1}]
decomposes as:

    J_{jk} = Σ_i J^(i)_{jk}  +  Σ_{i} W_cumsum_j(t_i) * dW^(i)_k

where J^(i) is the within-sub-interval iterated integral computed from the
dZ_i Fourier coefficients (Lévy area), and the second sum gives the
cross-interval contributions. Both terms are needed for correct strong
order 1.0 convergence on rejected steps.
"""
function _compute_II_from_S2(W_noise, m, dt)
    !hasproperty(W_noise, :S₂) && return nothing
    S₂ = W_noise.S₂
    n_sub = length(S₂)
    n_sub < 2 && return nothing

    T = Float64
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for idx in 1:n_sub
        dt_i = S₂.data[idx][1]  # sub-interval step size
        dWn = S₂.data[idx][2]   # dW_i from (dt_i, dW_i, dZ_i) tuple
        dZn = S₂.data[idx][3]   # dZ_i Fourier coefficients

        # Within-sub-interval iterated integral from Lévy area coefficients
        if dZn !== nothing && length(dZn) >= 2 * m && dt_i > 0
            coeffs_i = _unpack_dZ_to_coefficients(dZn, m, nothing)
            if coeffs_i !== nothing
                Wn_i = dWn / √dt_i
                A_i = levyarea(Wn_i, coeffs_i.n, MronRoe(), coeffs_i)
                # J^(i) = ½ dW_i dW_i' + dt_i * A_i  (Stratonovich within sub-interval)
                for k in 1:m
                    for j in 1:m
                        I[j, k] += 1 // 2 * dWn[j] * dWn[k] + dt_i * A_i[j, k]
                    end
                end
            else
                # No usable coefficients — add ½ dW_i dW_i' as best approximation
                for k in 1:m
                    for j in 1:m
                        I[j, k] += 1 // 2 * dWn[j] * dWn[k]
                    end
                end
            end
        else
            # No dZ available for this sub-interval — commutative approximation
            for k in 1:m
                for j in 1:m
                    I[j, k] += 1 // 2 * dWn[j] * dWn[k]
                end
            end
        end

        # Cross-interval contribution: W_cumsum_j * dW^(i)_k
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    return I
end

"""
Compute Stratonovich iterated integrals from saved fine-grid W and Z values in a
NoiseWrapper or NoiseProcess source. Uses Z grid (Fourier coefficients) for
within-sub-interval Lévy area when available, falling back to commutative
approximation (½ dW dW') for sub-intervals without Z data.
"""
function _compute_II_from_grid(W_noise, m, dt)
    source = _get_noise_grid_source(W_noise)
    source === nothing && return nothing

    t_grid = source.t
    W_grid = source.W
    t = W_noise.curt
    t_end = t + dt

    i_start = searchsortedfirst(t_grid, t)
    i_end = searchsortedlast(t_grid, t_end)

    # Need at least 2 sub-steps
    n_sub = i_end - i_start
    n_sub < 2 && return nothing

    # Check if Z grid is available for Lévy area computation
    Z_grid = _get_Z_grid(source)

    T = eltype(eltype(W_grid))
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for n in (i_start + 1):i_end
        dWn = W_grid[n] .- W_grid[n - 1]
        dt_n = t_grid[n] - t_grid[n - 1]

        # Within-sub-interval iterated integral
        if Z_grid !== nothing && dt_n > 0
            dZn = Z_grid[n] .- Z_grid[n - 1]
            if length(dZn) >= 2 * m
                coeffs_n = _unpack_dZ_to_coefficients(dZn, m, nothing)
                if coeffs_n !== nothing
                    Wn_n = dWn / √dt_n
                    A_n = levyarea(Wn_n, coeffs_n.n, MronRoe(), coeffs_n)
                    for k in 1:m
                        for j in 1:m
                            I[j, k] += 1 // 2 * dWn[j] * dWn[k] + dt_n * A_n[j, k]
                        end
                    end
                else
                    for k in 1:m
                        for j in 1:m
                            I[j, k] += 1 // 2 * dWn[j] * dWn[k]
                        end
                    end
                end
            else
                for k in 1:m
                    for j in 1:m
                        I[j, k] += 1 // 2 * dWn[j] * dWn[k]
                    end
                end
            end
        else
            # No Z data — commutative approximation
            for k in 1:m
                for j in 1:m
                    I[j, k] += 1 // 2 * dWn[j] * dWn[k]
                end
            end
        end

        # Cross-interval contribution
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    return I
end

"""
Extract Z grid from the noise source, or nothing if Z is not stored.
"""
function _get_Z_grid(source)
    hasproperty(source, :Z) || return nothing
    Z = source.Z
    (Z === nothing || isempty(Z)) && return nothing
    return Z
end

function _get_noise_grid_source(W)
    if hasproperty(W, :source)
        return _get_noise_grid_source(W.source)
    elseif W isa DiffEqNoiseProcess.NoiseGrid || W isa DiffEqNoiseProcess.NoiseProcess
        return W
    else
        return nothing
    end
end

@muladd function perform_step!(integrator, cache::RKMilGeneralConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    Jalg = cache.Jalg
    dW = W.dW

    # Try coefficient-based computation first (handles adaptivity + consistency)
    J = _compute_iterated_I(dt, dW, W.dZ, W, integrator.alg)
    if J === nothing
        # Fallback: legacy LevyArea path (scalar/diagonal/commutative noise)
        J = get_iterated_I(
            dt, dW, W.dZ, Jalg, integrator.alg.p, integrator.alg.c, alg_order(integrator.alg)
        )
    end

    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            J = J .- 1 // 2 .* abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    du1 = integrator.f(uprev, p, t)
    L = integrator.f.g(uprev, p, t)
    mil_correction = zero(u)
    ggprime_norm = zero(eltype(u))

    if dW isa Number || is_diagonal_noise(integrator.sol.prob)
        K = @.. uprev + dt * du1
        utilde = (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        mil_correction = ggprime .* J
        u = K + L .* dW + mil_correction
    else
        for i in 1:length(dW)
            K = uprev + dt * du1 + integrator.sqdt * @view(L[:, i])
            gtmp = integrator.f.g(K, p, t)
            ggprime = @.. (gtmp - L) / integrator.sqdt
            ggprime_norm = zero(eltype(u))
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(ggprime, t)
            end
            mil_correction += ggprime * @view(J[i, :])
        end
        if integrator.opts.adaptive
            K = @.. uprev + dt * du1
            u = K + L * dW + mil_correction
        else
            u = uprev + dt * du1 + L * dW + mil_correction
        end
    end

    if integrator.opts.adaptive
        du2 = integrator.f(K, p, t + dt)
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            tmp = dt * (du2 - du1) / 2
            En = W.dW .^ 3 .* ((du2 - L) / (integrator.sqdt)) .^ 2 / 6
        else
            En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
            tmp = integrator.opts.internalnorm((@.. dt * (du2 - du1) / 2), t)
        end
        tmp = calculate_residuals(
            tmp, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKMilGeneralCache)
    (; du1, du2, K, tmp, ggprime, L, mil_correction) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    dW = W.dW
    sqdt = integrator.sqdt
    Jalg = cache.Jalg

    # Try coefficient-based computation first
    J_coeffs = _compute_iterated_I(dt, dW, W.dZ, W, integrator.alg)
    if J_coeffs !== nothing
        Jalg.J .= J_coeffs
    else
        get_iterated_I!(
            dt, dW, W.dZ, Jalg, integrator.alg.p, integrator.alg.c, alg_order(integrator.alg)
        )
    end
    J = Jalg.J

    integrator.f(du1, uprev, p, t)
    integrator.f.g(L, uprev, p, t)
    @.. mil_correction = zero(eltype(u))
    ggprime_norm = zero(eltype(ggprime))

    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            @.. J -= 1 // 2 * abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    if dW isa Number || is_diagonal_noise(integrator.sol.prob)
        @.. K = uprev + dt * du1
        @.. du2 = zero(eltype(u))
        tmp .= (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) .+ integrator.sqdt .* L
        integrator.f.g(du2, tmp, p, t)
        @.. ggprime = (du2 - L) / sqdt
        ggprime_norm = integrator.opts.internalnorm(ggprime, t)
        @.. u = K + L * dW + ggprime * J
    else
        for i in 1:length(dW)
            @.. K = uprev + dt * du1 + sqdt * @view(L[:, i])
            integrator.f.g(ggprime, K, p, t)
            @.. ggprime = (ggprime - L) / sqdt
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(ggprime, t)
            end
            mul!(tmp, ggprime, @view(J[i, :]))
            @.. mil_correction += tmp
        end
        mul!(tmp, L, dW)
        if integrator.opts.adaptive
            @.. K = uprev + dt * du1
            @.. u = K + tmp + mil_correction
        else
            @.. u = uprev + dt * du1 + tmp + mil_correction
        end
    end

    if integrator.opts.adaptive
        En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
        integrator.f(du2, K, p, t + dt)
        @.. tmp = integrator.opts.internalnorm(
            integrator.opts.delta * dt * (du2 - du1) /
                2, t
        ) + En

        calculate_residuals!(
            tmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
    return nothing
end

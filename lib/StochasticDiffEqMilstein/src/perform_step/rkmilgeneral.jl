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

1. **NoiseWrapper/NoiseGrid sub-grid**: If the noise process has fine-grid W values
   (from a prior solve), compute iterated integrals from sub-grid dW increments.
   This ensures consistency between fine and coarse solutions for convergence testing.

2. **Rejection retry**: If dt < cached original dt, use sub-interval computation
   from the cached Fourier coefficients.

3. **Fresh step with dZ coefficients**: Unpack dZ into MronRoe Fourier coefficients,
   compute levyarea deterministically. Cache for potential rejection.

4. **Fallback**: Return nothing → legacy LevyArea path.
"""
function _compute_iterated_I(dt, dW, dZ, W_noise, cache, alg)
    m = length(dW)

    # Check if we have usable dZ for coefficient-based computation
    if dZ === nothing || length(dZ) < 2 * m
        return nothing
    end

    # Strategy 1: NoiseWrapper/NoiseGrid sub-grid (convergence testing)
    J_subgrid = _compute_II_from_grid(W_noise, m, dt)
    if J_subgrid !== nothing
        return J_subgrid
    end

    # Strategy 2: Rejection retry — use cached original coefficients
    if cache._dt_orig[] > 0 && dt < cache._dt_orig[] - eps(cache._dt_orig[])
        dZ_orig = cache._dZ_orig[]
        dW_orig = cache._dW_orig[]
        dt_orig = cache._dt_orig[]
        coeffs = _unpack_dZ_to_coefficients(dZ_orig, m, alg.p)
        coeffs === nothing && return nothing

        n_quad = max(64, coeffs.n)
        J = iterated_integrals_subinterval(
            dW_orig, dt_orig, coeffs, 0.0, dt;
            n_quadrature = n_quad, ito_correction = false
        )
        return J
    end

    # Strategy 3: Fresh step — unpack dZ into Fourier coefficients
    coeffs = _unpack_dZ_to_coefficients(dZ, m, alg.p)
    coeffs === nothing && return nothing

    # Cache for potential rejection retry
    copyto!(cache._dW_orig[], dW)
    copyto!(cache._dZ_orig[], dZ)
    cache._dt_orig[] = dt

    # Compute full-step Lévy area from coefficients (Stratonovich)
    Wn = dW / √dt
    A = levyarea(Wn, coeffs.n, MronRoe(), coeffs)
    J = 1 // 2 * dW .* dW' .+ dt .* A
    return J
end

"""
Compute Stratonovich iterated integrals from saved fine-grid W values in a
NoiseWrapper or NoiseProcess source. Returns nothing if no sub-grid is available.
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

    T = eltype(eltype(W_grid))
    I = zeros(T, m, m)
    W_cumsum = zeros(T, m)

    for n in (i_start + 1):i_end
        dWn = W_grid[n] .- W_grid[n - 1]
        for k in 1:m
            for j in 1:m
                I[j, k] += W_cumsum[j] * dWn[k]
            end
        end
        W_cumsum .+= dWn
    end

    # Ito → Stratonovich
    for j in 1:m
        I[j, j] += dt / 2
    end
    return I
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
    J = _compute_iterated_I(dt, dW, W.dZ, W, cache, integrator.alg)
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
    J_coeffs = _compute_iterated_I(dt, dW, W.dZ, W, cache, integrator.alg)
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

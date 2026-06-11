# Stage 1 is identity (Y_1 = u_n); inner ODE active only for i ≥ 2.

function initialize!(integrator, cache::MIS2Cache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst .+= cache.tmp
end

function initialize!(integrator, cache::MIS2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::MIS2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, offset, k_fast, Y, fS, tab) = cache
    (; α, β, γ, d, c, ctilde) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(d)

    @.. broadcast = false Y[1] = uprev
    f.f2(fS[1], uprev, p, t)
    integrator.stats.nf2 += 1

    for i in 2:s
        d_i = d[i]
        inv_di = 1 / d_i

        @.. broadcast = false v = uprev
        for j in 1:(i - 1)
            if !iszero(α[i, j])
                @.. broadcast = false v = v + α[i, j] * (Y[j] - uprev)
            end
        end

        @.. broadcast = false offset = zero(offset)
        for j in 1:(i - 1)
            if !iszero(γ[i, j])
                @.. broadcast = false offset = offset + (γ[i, j] * inv_di / dt) * (Y[j] - uprev)
            end
            if !iszero(β[i, j])
                @.. broadcast = false offset = offset + (β[i, j] * inv_di) * fS[j]
            end
        end

        M_inner = max(1, ceil(Int, m * d_i))
        h_inner = d_i * dt / M_inner
        slope = (c[i] - ctilde[i]) * inv_di
        for k_step in 1:M_inner
            τ = (k_step - 1) * h_inner
            t_a = t + ctilde[i] * dt + slope * τ
            t_b = t + ctilde[i] * dt + slope * (τ + h_inner / 2)
            f.f1(k_fast, v, p, t_a)
            @.. broadcast = false tmp = v + (h_inner / 2) * (offset + k_fast)
            f.f1(k_fast, tmp, p, t_b)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            @.. broadcast = false v = v + h_inner * (offset + k_fast)
        end

        @.. broadcast = false Y[i] = v
        f.f2(fS[i], Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1
    end

    @.. broadcast = false u = Y[s]

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = Y[s] - Y[s - 1]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(integrator, cache::MIS2ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    (; α, β, γ, d, c, ctilde) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(d)

    Y = Vector{typeof(uprev)}(undef, s)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)

    Y[1] = uprev
    fS[1] = f.f2(uprev, p, t)
    integrator.stats.nf2 += 1

    for i in 2:s
        d_i = d[i]
        inv_di = 1 / d_i

        v_state = uprev
        for j in 1:(i - 1)
            if !iszero(α[i, j])
                v_state = @.. broadcast = false v_state + α[i, j] * (Y[j] - uprev)
            end
        end

        offset = zero(fS[1])
        for j in 1:(i - 1)
            if !iszero(γ[i, j])
                offset = @.. broadcast = false offset + (γ[i, j] * inv_di / dt) * (Y[j] - uprev)
            end
            if !iszero(β[i, j])
                offset = @.. broadcast = false offset + (β[i, j] * inv_di) * fS[j]
            end
        end

        M_inner = max(1, ceil(Int, m * d_i))
        h_inner = d_i * dt / M_inner
        slope = (c[i] - ctilde[i]) * inv_di
        for k_step in 1:M_inner
            τ = (k_step - 1) * h_inner
            t_a = t + ctilde[i] * dt + slope * τ
            t_b = t + ctilde[i] * dt + slope * (τ + h_inner / 2)
            k_fast = f.f1(v_state, p, t_a)
            v_mid = @.. broadcast = false v_state + (h_inner / 2) * (offset + k_fast)
            k_fast = f.f1(v_mid, p, t_b)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            v_state = @.. broadcast = false v_state + h_inner * (offset + k_fast)
        end

        Y[i] = v_state
        fS[i] = f.f2(Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1
    end

    integrator.u = Y[s]

    if integrator.opts.adaptive
        utilde = @.. broadcast = false Y[s] - Y[s - 1]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

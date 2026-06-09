function initialize!(integrator, cache::MRIGARKERK22Cache)
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

function initialize!(integrator, cache::MRIGARKERK22ConstantCache)
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

function perform_step!(integrator, cache::MRIGARKERK22Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, f1eval, fS_yn, fS_Y2, Y2, tab) = cache
    (; Γ11, Γ22, γ11, γ21, γ22) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m

    f.f2(fS_yn, uprev, p, t)
    integrator.stats.nf2 += 1

    @.. broadcast = false v = uprev
    h_inner = dt / m
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        f.f1(f1eval, v, p, t_a)
        @.. broadcast = false tmp = v + (h_inner / 2) * (Γ11 * f1eval + γ11 * fS_yn)
        f.f1(f1eval, tmp, p, t_a + h_inner / 2)
        @.. broadcast = false v = v + h_inner * (Γ11 * f1eval + γ11 * fS_yn)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    @.. broadcast = false Y2 = v

    f.f2(fS_Y2, Y2, p, t + dt)
    integrator.stats.nf2 += 1

    if iszero(Γ22)
        @.. broadcast = false u = Y2 + dt * (γ21 * fS_yn + γ22 * fS_Y2)
    else
        @.. broadcast = false v = Y2
        for k_step in 1:m
            t_a = t + (k_step - 1) * h_inner
            f.f1(f1eval, v, p, t_a)
            @.. broadcast = false tmp = v + (h_inner / 2) *
                (Γ22 * f1eval + γ21 * fS_yn + γ22 * fS_Y2)
            f.f1(f1eval, tmp, p, t_a + h_inner / 2)
            @.. broadcast = false v = v + h_inner *
                (Γ22 * f1eval + γ21 * fS_yn + γ22 * fS_Y2)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        end
        @.. broadcast = false u = v
    end

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = u - Y2
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(
        integrator, cache::MRIGARKERK22ConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; Γ11, Γ22, γ11, γ21, γ22) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m

    fS_yn = f.f2(uprev, p, t)
    integrator.stats.nf2 += 1
    h_inner = dt / m

    v = uprev
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        k1 = f.f1(v, p, t_a)
        v_mid = @.. broadcast = false v + (h_inner / 2) * (Γ11 * k1 + γ11 * fS_yn)
        k2 = f.f1(v_mid, p, t_a + h_inner / 2)
        v = @.. broadcast = false v + h_inner * (Γ11 * k2 + γ11 * fS_yn)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    Y2 = v

    fS_Y2 = f.f2(Y2, p, t + dt)
    integrator.stats.nf2 += 1

    if iszero(Γ22)
        integrator.u = @.. broadcast = false Y2 + dt * (γ21 * fS_yn + γ22 * fS_Y2)
    else
        v = Y2
        for k_step in 1:m
            t_a = t + (k_step - 1) * h_inner
            k1 = f.f1(v, p, t_a)
            v_mid = @.. broadcast = false v + (h_inner / 2) *
                (Γ22 * k1 + γ21 * fS_yn + γ22 * fS_Y2)
            k2 = f.f1(v_mid, p, t_a + h_inner / 2)
            v = @.. broadcast = false v + h_inner *
                (Γ22 * k2 + γ21 * fS_yn + γ22 * fS_Y2)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        end
        integrator.u = v
    end

    if integrator.opts.adaptive
        utilde = @.. broadcast = false integrator.u - Y2
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

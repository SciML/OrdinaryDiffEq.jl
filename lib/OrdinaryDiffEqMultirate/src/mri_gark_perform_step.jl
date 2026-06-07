# MRI-GARK-ERK22a — Sandu, "A class of Multirate Infinitesimal GARK methods",
# SIAM J. Numer. Anal. 57 (2019), 2300–2327, equation (2.8).
#
# For a SplitODE du/dt = f1(u, t) + f2(u, t)  (f1 fast, f2 slow):
#
#   Stage 1:  v_1(0) = u_n;   dv_1/dτ = (1/2) f1(v_1) + (1/2) f2(u_n)
#             integrate τ ∈ [0, dt], Y_2 = v_1(dt)
#   Stage 2:  v_2(0) = Y_2;   dv_2/dτ = (1/2) f1(v_2) - (1/2) f2(u_n) + f2(Y_2)
#             integrate τ ∈ [0, dt], u_{n+1} = v_2(dt)
#
# Inner micro-ODE integrated with explicit-midpoint (RK2) micro-steps to keep
# the inner contribution from capping the design order 2.

function initialize!(integrator, cache::MRIGARKERK22aCache)
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

function initialize!(integrator, cache::MRIGARKERK22aConstantCache)
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

function perform_step!(integrator, cache::MRIGARKERK22aCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, f1eval, fS_yn, fS_Y2, Y2) = cache
    alg = unwrap_alg(integrator, false)
    m = alg.m

    f.f2(fS_yn, uprev, p, t)
    integrator.stats.nf2 += 1

    @.. broadcast = false v = uprev
    h_inner = dt / m
    offset1 = fS_yn
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        f.f1(f1eval, v, p, t_a)
        @.. broadcast = false tmp = v + (h_inner / 2) * ((1 // 2) * f1eval + (1 // 2) * offset1)
        f.f1(f1eval, tmp, p, t_a + h_inner / 2)
        @.. broadcast = false v = v + h_inner * ((1 // 2) * f1eval + (1 // 2) * offset1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    @.. broadcast = false Y2 = v

    f.f2(fS_Y2, Y2, p, t + dt)
    integrator.stats.nf2 += 1

    @.. broadcast = false v = Y2
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        f.f1(f1eval, v, p, t_a)
        @.. broadcast = false tmp = v + (h_inner / 2) *
            ((1 // 2) * f1eval - (1 // 2) * fS_yn + fS_Y2)
        f.f1(f1eval, tmp, p, t_a + h_inner / 2)
        @.. broadcast = false v = v + h_inner *
            ((1 // 2) * f1eval - (1 // 2) * fS_yn + fS_Y2)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    @.. broadcast = false u = v

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = v - Y2
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(integrator, cache::MRIGARKERK22aConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    m = alg.m

    fS_yn = f.f2(uprev, p, t)
    integrator.stats.nf2 += 1

    h_inner = dt / m

    v = uprev
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        k1 = f.f1(v, p, t_a)
        v_mid = @.. broadcast = false v + (h_inner / 2) * ((1 // 2) * k1 + (1 // 2) * fS_yn)
        k2 = f.f1(v_mid, p, t_a + h_inner / 2)
        v = @.. broadcast = false v + h_inner * ((1 // 2) * k2 + (1 // 2) * fS_yn)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    Y2 = v

    fS_Y2 = f.f2(Y2, p, t + dt)
    integrator.stats.nf2 += 1

    v = Y2
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        k1 = f.f1(v, p, t_a)
        v_mid = @.. broadcast = false v + (h_inner / 2) *
            ((1 // 2) * k1 - (1 // 2) * fS_yn + fS_Y2)
        k2 = f.f1(v_mid, p, t_a + h_inner / 2)
        v = @.. broadcast = false v + h_inner *
            ((1 // 2) * k2 - (1 // 2) * fS_yn + fS_Y2)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    integrator.u = v

    if integrator.opts.adaptive
        utilde = @.. broadcast = false v - Y2
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

# MRI-GARK-ERK22b — Sandu 2019, eq. (2.7) family with c₂ = 1 (explicit trapezoidal).
#
# For c₂ = 1 the second-stage fast coefficient Γ₂₂ = 1 - c₂ = 0, so stage 2
# carries no fast dynamics; it is a pure slow correction:
#
#   Stage 1: v(0) = u_n;  dv/dτ = f1(v) + f2(u_n)  τ ∈ [0,dt]  → Y_2 = v(dt)
#   Stage 2: u_{n+1} = Y_2 + dt·(−½ f2(u_n) + ½ f2(Y_2))

function initialize!(integrator, cache::MRIGARKERK22bCache)
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

function initialize!(integrator, cache::MRIGARKERK22bConstantCache)
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

function perform_step!(integrator, cache::MRIGARKERK22bCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, f1eval, fS_yn, fS_Y2, Y2) = cache
    alg = unwrap_alg(integrator, false)
    m = alg.m

    f.f2(fS_yn, uprev, p, t)
    integrator.stats.nf2 += 1

    @.. broadcast = false v = uprev
    h_inner = dt / m
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        f.f1(f1eval, v, p, t_a)
        @.. broadcast = false tmp = v + (h_inner / 2) * (f1eval + fS_yn)
        f.f1(f1eval, tmp, p, t_a + h_inner / 2)
        @.. broadcast = false v = v + h_inner * (f1eval + fS_yn)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    @.. broadcast = false Y2 = v

    f.f2(fS_Y2, Y2, p, t + dt)
    integrator.stats.nf2 += 1

    @.. broadcast = false u = Y2 + dt * (-(1 // 2) * fS_yn + (1 // 2) * fS_Y2)

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
        integrator, cache::MRIGARKERK22bConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    m = alg.m

    fS_yn = f.f2(uprev, p, t)
    integrator.stats.nf2 += 1
    h_inner = dt / m

    v = uprev
    for k_step in 1:m
        t_a = t + (k_step - 1) * h_inner
        k1 = f.f1(v, p, t_a)
        v_mid = @.. broadcast = false v + (h_inner / 2) * (k1 + fS_yn)
        k2 = f.f1(v_mid, p, t_a + h_inner / 2)
        v = @.. broadcast = false v + h_inner * (k2 + fS_yn)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    end
    Y2 = v

    fS_Y2 = f.f2(Y2, p, t + dt)
    integrator.stats.nf2 += 1

    integrator.u = @.. broadcast = false Y2 + dt * (-(1 // 2) * fS_yn + (1 // 2) * fS_Y2)

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

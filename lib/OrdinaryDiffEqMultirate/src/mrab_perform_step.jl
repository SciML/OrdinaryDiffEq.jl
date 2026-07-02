# ── MRAB: Adams–Bashforth coefficients ────────────────────────────────────────
#
# Standard explicit AB-k weights `β` ordered so that
#   u_{n+1} = u_n + h * Σ_{i=1}^k β[i] * F_{n+1-i}
# (β[1] applies to the most recent rate sample).

@inline function _ab_betas(k::Int)
    if k == 1
        return (1.0,)
    elseif k == 2
        return (3 / 2, -1 / 2)
    elseif k == 3
        return (23 / 12, -16 / 12, 5 / 12)
    elseif k == 4
        return (55 / 24, -59 / 24, 37 / 24, -9 / 24)
    elseif k == 5
        return (1901 / 720, -2774 / 720, 2616 / 720, -1274 / 720, 251 / 720)
    else
        throw(ArgumentError("MRAB: unsupported Adams order k=$k (supported: 1..5)"))
    end
end

# ── MRAB initialize! ──────────────────────────────────────────────────────────

function initialize!(integrator, cache::MRABCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    (; tmp) = cache.tmp_cache
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)  # f1
    integrator.stats.nf2 += 1  # f2
    return integrator.fsalfirst .+= tmp
end

function initialize!(integrator, cache::MRABConstantCache)
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

# ── MRAB perform_step! (in-place, MutableCache) ───────────────────────────────
# Substep ℓ bootstraps with AB-min(ℓ, k) before the history reaches k samples.

function perform_step!(integrator, cache::MRABCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k_slow, k_fast, F_history) = cache
    (; tmp, atmp) = cache.tmp_cache
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = _ab_betas(k)

    @.. broadcast = false u = uprev
    f.f2(k_slow, u, p, t)
    integrator.stats.nf2 += 1

    n_hist = 0
    for ℓ in 1:m
        t_fast = t + (ℓ - 1) * h
        f.f1(k_fast, u, p, t_fast)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        @.. broadcast = false F_history[min(n_hist + 1, k)] = k_fast + k_slow
        for j in min(n_hist + 1, k):-1:2
            F_history[j], F_history[j - 1] = F_history[j - 1], F_history[j]
        end
        n_hist = min(n_hist + 1, k)

        ks = min(ℓ, k)
        βs_use = _ab_betas(ks)
        for i in 1:ks
            β = βs_use[i]
            Fi = F_history[i]
            @.. broadcast = false u = u + h * β * Fi
        end
    end

    # Error estimate (AB-k − AB-(k-1)) folds into one linear combination over F_history.
    return if integrator.opts.adaptive && k >= 2
        βs_low = _ab_betas(k - 1)
        @.. broadcast = false tmp = h * (βs[1] - βs_low[1]) * F_history[1]
        for i in 2:(k - 1)
            @.. broadcast = false tmp = tmp + h * (βs[i] - βs_low[i]) * F_history[i]
        end
        @.. broadcast = false tmp = tmp + h * βs[k] * F_history[k]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

# ── MRAB perform_step! (out-of-place, ConstantCache) ──────────────────────────

@muladd function perform_step!(integrator, cache::MRABConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = _ab_betas(k)

    u_cur = uprev
    k_slow = f.f2(u_cur, p, t)
    integrator.stats.nf2 += 1

    F_history = Vector{typeof(k_slow)}(undef, k)
    n_hist = 0

    for ℓ in 1:m
        t_fast = t + (ℓ - 1) * h
        k_fast = f.f1(u_cur, p, t_fast)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        F = k_fast + k_slow
        n_hist = min(n_hist + 1, k)
        for j in n_hist:-1:2
            F_history[j] = F_history[j - 1]
        end
        F_history[1] = F

        ks = min(ℓ, k)
        βs_use = _ab_betas(ks)
        Δu = zero(u_cur)
        for i in 1:ks
            Δu = @.. broadcast = false Δu + h * βs_use[i] * F_history[i]
        end
        u_cur = @.. broadcast = false u_cur + Δu
    end

    integrator.u = u_cur

    if integrator.opts.adaptive && k >= 2
        βs_low = _ab_betas(k - 1)
        utilde = @.. broadcast = false h * (βs[1] - βs_low[1]) * F_history[1]
        for i in 2:(k - 1)
            utilde = @.. broadcast = false utilde + h * (βs[i] - βs_low[i]) * F_history[i]
        end
        utilde = @.. broadcast = false utilde + h * βs[k] * F_history[k]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

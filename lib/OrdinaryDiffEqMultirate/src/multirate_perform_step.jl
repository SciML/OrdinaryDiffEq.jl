# ── MREEF: step-count sequence ────────────────────────────────────────────────

@inline function _mreef_sequence(seq::Symbol, order::Int)
    if seq === :harmonic
        return ntuple(j -> j, order)
    elseif seq === :romberg
        return ntuple(j -> 1 << (j - 1), order)
    else
        throw(ArgumentError("MREEF: unknown sequence `$seq`, choose :harmonic or :romberg"))
    end
end

# ── MREEF initialize! ─────────────────────────────────────────────────────────

function initialize!(integrator, cache::MREEFCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)  # f1
    integrator.stats.nf2 += 1  # f2
    return integrator.fsalfirst .+= cache.tmp
end

function initialize!(integrator, cache::MREEFConstantCache)
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

# ── MREEF perform_step! (in-place, MutableCache) ──────────────────────────────
#
# Base multirate Euler with nj macro intervals, m fast substeps each:
#   1. k_slow = f.f2(u, p, t_mac)  — frozen slow rate for the macro interval
#   2. m fast substeps: u += h_fast*(k_slow + f.f1(u, p, t_fast))
# f1 = fast/stiff (large eigenvalues), f2 = slow/non-stiff (SciML convention).
# Then apply Aitken–Neville Richardson extrapolation over T[1..order].

function perform_step!(integrator, cache::MREEFCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, k_slow, k_fast, T) = cache
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _mreef_sequence(alg.seq, order)

    # Fill first tableau column: T[j] = base method with ns[j] macro intervals
    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        @.. broadcast = false T[j] = uprev

        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac

            # Slow evaluation (f2): frozen for all m fast substeps
            f.f2(k_slow, T[j], p, t_mac)
            integrator.stats.nf2 += 1

            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                f.f1(k_fast, T[j], p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast = false T[j] = T[j] + h_fast * k_slow + h_fast * k_fast
            end
        end
    end

    # Aitken–Neville Richardson extrapolation (in-place, reverse-row order)
    # Formula: T[j] <- T[j] + (T[j] - T[j-1]) / (ns[j]/ns[j-k] - 1)
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
            @.. broadcast = false tmp = (T[j] - T[j - 1]) / (ratio - 1)
            @.. broadcast = false T[j] = T[j] + tmp
        end
    end

    @.. broadcast = false u = T[order]

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = T[order] - T[order - 1]
        calculate_residuals!(
            atmp,
            tmp,
            uprev,
            u,
            integrator.opts.abstol,
            integrator.opts.reltol,
            integrator.opts.internalnorm,
            t,
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

# ── MREEF perform_step! (out-of-place, ConstantCache) ─────────────────────────

@muladd function perform_step!(integrator, cache::MREEFConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _mreef_sequence(alg.seq, order)
    T = cache.T

    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        u_cur = uprev
        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac
            k_slow = f.f2(u_cur, p, t_mac)
            integrator.stats.nf2 += 1
            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                k_fast = f.f1(u_cur, p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                u_cur = @.. broadcast = false u_cur + h_fast * k_slow + h_fast * k_fast
            end
        end
        T[j] = u_cur
    end

    # Aitken–Neville Richardson extrapolation
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
            T[j] = @.. broadcast = false T[j] + (T[j] - T[j - 1]) / (ratio - 1)
        end
    end

    integrator.u = T[order]

    if integrator.opts.adaptive
        utilde = @.. broadcast = false T[order] - T[order - 1]
        atmp = calculate_residuals(
            utilde,
            uprev,
            integrator.u,
            integrator.opts.abstol,
            integrator.opts.reltol,
            integrator.opts.internalnorm,
            t,
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

function initialize!(integrator, cache::MRABCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)  # f1
    integrator.stats.nf2 += 1  # f2
    return integrator.fsalfirst .+= cache.tmp
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
    (; tmp, atmp, k_slow, k_fast, F_history, tab) = cache
    β = tab.β
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = β[k]

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
        βs_use = β[ks]
        for i in 1:ks
            βi = βs_use[i]
            Fi = F_history[i]
            @.. broadcast = false u = u + h * βi * Fi
        end
    end

    # Error estimate (AB-k − AB-(k-1)) folds into one linear combination over F_history.
    return if integrator.opts.adaptive && k >= 2
        βs_low = β[k - 1]
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
    β = cache.tab.β
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = β[k]

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
        βs_use = β[ks]
        Δu = zero(u_cur)
        for i in 1:ks
            Δu = @.. broadcast = false Δu + h * βs_use[i] * F_history[i]
        end
        u_cur = @.. broadcast = false u_cur + Δu
    end

    integrator.u = u_cur

    if integrator.opts.adaptive && k >= 2
        βs_low = β[k - 1]
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

function initialize!(integrator, cache::MRIGARKCache)
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

function initialize!(integrator, cache::MRIGARKConstantCache)
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

function perform_step!(integrator, cache::MRIGARKCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, f1eval, slow_f, Y, fS, tab) = cache
    (; Γ, γ, c) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Γ)

    @.. broadcast = false Y[1] = uprev
    h_inner = dt / m

    for i in 1:s
        f.f2(fS[i], Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1

        @.. broadcast = false slow_f = γ[i, 1] * fS[1]
        for j in 2:i
            if !iszero(γ[i, j])
                @.. broadcast = false slow_f = slow_f + γ[i, j] * fS[j]
            end
        end

        if iszero(Γ[i])
            if i < s
                @.. broadcast = false Y[i + 1] = Y[i] + dt * slow_f
            else
                @.. broadcast = false u = Y[i] + dt * slow_f
            end
        else
            @.. broadcast = false v = Y[i]
            for k_step in 1:m
                t_a = t + (k_step - 1) * h_inner
                f.f1(f1eval, v, p, t_a)
                @.. broadcast = false tmp = v + (h_inner / 2) * (Γ[i] * f1eval + slow_f)
                f.f1(f1eval, tmp, p, t_a + h_inner / 2)
                @.. broadcast = false v = v + h_inner * (Γ[i] * f1eval + slow_f)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            end
            if i < s
                @.. broadcast = false Y[i + 1] = v
            else
                @.. broadcast = false u = v
            end
        end
    end

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = u - Y[s]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(
        integrator, cache::MRIGARKConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; Γ, γ, c) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Γ)

    h_inner = dt / m
    Y = Vector{typeof(uprev)}(undef, s)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)
    Y[1] = uprev

    for i in 1:s
        fS[i] = f.f2(Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1

        slow_f = γ[i, 1] * fS[1]
        for j in 2:i
            if !iszero(γ[i, j])
                slow_f = @.. broadcast = false slow_f + γ[i, j] * fS[j]
            end
        end

        if iszero(Γ[i])
            result = @.. broadcast = false Y[i] + dt * slow_f
        else
            v = Y[i]
            for k_step in 1:m
                t_a = t + (k_step - 1) * h_inner
                k1 = f.f1(v, p, t_a)
                v_mid = @.. broadcast = false v + (h_inner / 2) * (Γ[i] * k1 + slow_f)
                k2 = f.f1(v_mid, p, t_a + h_inner / 2)
                v = @.. broadcast = false v + h_inner * (Γ[i] * k2 + slow_f)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            end
            result = v
        end

        if i < s
            Y[i + 1] = result
        else
            integrator.u = result
        end
    end

    if integrator.opts.adaptive
        utilde = @.. broadcast = false integrator.u - Y[s]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

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

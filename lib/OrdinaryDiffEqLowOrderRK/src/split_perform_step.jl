function initialize!(integrator, cache::SplitEulerConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::SplitEulerConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    u = @.. broadcast = false uprev + dt * integrator.fsalfirst
    integrator.fsallast = f.f1(u, p, t + dt) + f.f2(u, p, t + dt)  # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

get_fsalfirstlast(cache::SplitEulerCache, u) = (cache.fsalfirst, cache.k)
function initialize!(integrator, cache::SplitEulerCache)
    integrator.kshortsize = 2
    (; k, fsalfirst) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst .+= cache.tmp
end

@muladd function perform_step!(integrator, cache::SplitEulerCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    @.. broadcast = false u = uprev + dt * integrator.fsalfirst
    f.f1(integrator.fsallast, u, p, t + dt) # For the interpolation, needs k at the updated point
    f.f2(cache.tmp, u, p, t + dt) # For the interpolation, needs k at the updated point
    integrator.stats.nf2 += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast .+= cache.tmp
end

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
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
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
#   1. k_slow = f.f1(u, p, t_mac)  — frozen slow rate for the macro interval
#   2. m fast substeps: u += h_fast*(k_slow + f.f2(u, p, t_fast))
# Then apply Aitken–Neville Richardson extrapolation over T[1..order].

function perform_step!(integrator, cache::MREEFCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, k_slow, k_fast, T) = cache
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _mreef_sequence(alg.seq, order)

    # Fill first tableau column: T[j] = base method with ns[j] macro intervals
    for j = 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        @.. broadcast=false T[j] = uprev

        for i_mac = 1:nj
            t_mac = t + (i_mac - 1) * h_mac

            # Slow evaluation: frozen for all m fast substeps
            f.f1(k_slow, T[j], p, t_mac)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

            for i_fast = 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                f.f2(k_fast, T[j], p, t_fast)
                integrator.stats.nf2 += 1
                @.. broadcast=false T[j] = T[j] + h_fast * k_slow + h_fast * k_fast
            end
        end
    end

    # Aitken–Neville Richardson extrapolation (in-place, reverse-row order)
    # Formula: T[j] <- T[j] + (T[j] - T[j-1]) / (ns[j]/ns[j-k] - 1)
    for k = 1:(order-1)
        for j = order:-1:(k+1)
            ratio = ns[j] / ns[j-k]
            @.. broadcast=false tmp = (T[j] - T[j-1]) / (ratio - 1)
            @.. broadcast=false T[j] = T[j] + tmp
        end
    end

    @.. broadcast=false u = T[order]

    if integrator.opts.adaptive
        @.. broadcast=false tmp = T[order] - T[order-1]
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
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

# ── MREEF perform_step! (out-of-place, ConstantCache) ─────────────────────────

@muladd function perform_step!(integrator, cache::MREEFConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _mreef_sequence(alg.seq, order)

    T = Vector{typeof(uprev)}(undef, order)

    for j = 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        u_cur = uprev
        for i_mac = 1:nj
            t_mac = t + (i_mac - 1) * h_mac
            k_slow = f.f1(u_cur, p, t_mac)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            for i_fast = 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                k_fast = f.f2(u_cur, p, t_fast)
                integrator.stats.nf2 += 1
                u_cur = @.. broadcast=false u_cur + h_fast * k_slow + h_fast * k_fast
            end
        end
        T[j] = u_cur
    end

    # Aitken–Neville Richardson extrapolation
    for k = 1:(order-1)
        for j = order:-1:(k+1)
            ratio = ns[j] / ns[j-k]
            T[j] = @.. broadcast=false T[j] + (T[j] - T[j-1]) / (ratio - 1)
        end
    end

    u = T[order]

    if integrator.opts.adaptive
        utilde = @.. broadcast=false T[order] - T[order-1]
        atmp = calculate_residuals(
            utilde,
            uprev,
            u,
            integrator.opts.abstol,
            integrator.opts.reltol,
            integrator.opts.internalnorm,
            t,
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

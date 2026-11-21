# Newmark-β is generalized-α with αm = αf = 0, and both schemes reach the inner solver
# through the same discretization cache, so the in-place step is shared.

# The step re-derives aₙ from f at the accepted state while aₙ₊₁ solves the
# interpolated-state balance M[(1 - αm)aₙ₊₁ + αm aₙ] = f(xₙ₊αf), which stretches
# the solved acceleration difference: aₙ₊₁ - aₙ = κ dt u''' + O(dt²) with
# κ = (1 - αf)/(1 - αm). Newmark-β has κ = 1.
zx_kappa(αm, αf) = (1 - αf) / (1 - αm)

# Zienkiewicz & Xie (1991) Eq. 21 generalized to that discretization: comparing
# the update formulas against the Taylor expansion of the exact solution, the
# leading local errors of a step are
#   position: (βκ - 1/6) dt³ u''' = (β - 1/(6κ)) dt² (aₙ₊₁ - aₙ)
#   velocity: (γκ - 1/2) dt² u''' = (γ - 1/(2κ)) dt  (aₙ₊₁ - aₙ)
# For Newmark-β the position line is Eq. 21 verbatim.
zx_position_coefficient(β, κ) = β - 1 / (6κ)
zx_velocity_coefficient(γ, κ) = γ - 1 / (2κ)

# Each channel reports an error arbitrarily smaller than the true one near its
# vanishing locus: β = 1/(6κ) for position, γ = 1/(2κ) for velocity. The two
# channels are not interchangeable. A healthy velocity channel does bound the
# position error, because it watches the lower-order dt² term and the position
# error is dt³ (measured on the harmonic oscillator: β = 1/6, γ = 0.7 — zero
# position coefficient — holds Q = global error / (nsteps * tol) ≈ 1.0 over
# tol = 1e-6..1e-10). A healthy position channel does NOT bound the velocity error:
# on the family γ = 1/(2κ) both leading errors are O(dt³), but the velocity
# term carries an extra factor of frequency (dt³ω³ against the position dt³ω²
# on u'' = -ω²u), so the position-equilibrated dt lets the velocity error pass
# the tolerance by a factor that grows linearly in ω — measured Qv ≈ 2.1/13.5/119
# at ω = 1/10/100 (tol = 1e-6) and full 2.0 peak-to-peak velocity error with a
# Success code at ω = 1000, before the clamp below existed.
#
# The velocity channel therefore never switches off. Inside the band
# |γ - 1/(2κ)| < ZX_VELOCITY_BAND the measured coefficient is exactly the
# quantity the band declares untrustworthy, so the channel is evaluated as if
# γ sat on the band edge: Cv_eff = max(|γ - 1/(2κ)|, ZX_VELOCITY_BAND). The
# in-band estimate Cv_eff·dt·(aₙ₊₁ - aₙ) ≈ (κ/100)·dt²·u''' is a surrogate one
# order below the true on-family velocity error (γκ·c₂ - 1/6)·dt³·u'''' — the
# next-order term is not measurable from a single step's Δa — hence
# conservative for resolved motion and continuous in γ across the band edge.
# The cost is that dt scales as tol^(1/2) instead of tol^(1/3) wherever the
# surrogate binds; measured step counts at ω = 1 on the harmonic oscillator:
# 1.05x/1.9x/4.2x/9.1x the position-only counts at tol = 1e-4/-6/-8/-10
# (43/186/838/3863 steps before, 45/361/3481/35023 after), and with the
# surrogate binding the omega sweep ω = 1..1000, tol = 1e-4..1e-10 measures
# Q ≤ 2.5 everywhere (worst measured ≈ 2.46 near ω = 30, tol = 1e-4).
#
# Both bands together (position channel vanishing AND γ in the velocity band)
# are still rejected: the estimate would consist solely of the surrogate,
# which is a proxy for step control, not a measurement of either leading
# error term. The position width 1/24 predates the surrogate (position-channel
# quality degraded like ~0.19/|Cp| when nothing else bound the step) and is
# kept as a conservative gate for that rejection. Baseline quality on the
# harmonic oscillator over reltol = abstol in [1e-12, 1e-4] is Q ≈ 0.2-3
# for every healthy configuration measured. Q is problem-dependent: the
# near-separatrix pendulum (u0 = 2.5 rad, T = 10) measures Q ≈ 34-39 for the
# velocity-channel configuration β = 1/6, γ = 0.7 at the same tolerances,
# while a small-angle pendulum (u0 = 0.5 rad) stays at Q ≈ 0.6-1.4.
const ZX_POSITION_BAND = 1 / 24
const ZX_VELOCITY_BAND = 1 / 100

function zx_band_check(β, γ, αm, αf)
    κ = zx_kappa(αm, αf)
    Cp = zx_position_coefficient(β, κ)
    Cv = zx_velocity_coefficient(γ, κ)
    if abs(Cp) < ZX_POSITION_BAND && abs(Cv) < ZX_VELOCITY_BAND
        throw(
            ArgumentError(
                "The Zienkiewicz-Xie error estimate is blind for these " *
                    "parameters: its position coefficient β - 1/(6κ) = $Cp and " *
                    "velocity coefficient γ - 1/(2κ) = $Cv, κ = (1 - αf)/(1 - αm), " *
                    "are both too close to their vanishing loci, so the estimated " *
                    "error can be arbitrarily small next to the actual error and " *
                    "the solver would return inaccurate results with a Success " *
                    "code. Move β at least $ZX_POSITION_BAND away from $(1 / (6κ)) " *
                    "or γ at least $ZX_VELOCITY_BAND away from $(1 / (2κ)), or " *
                    "pass adaptive = false."
            )
        )
    end
    return nothing
end

function initialize!(
        integrator, cache::Union{NewmarkBetaCache, GeneralizedAlphaCache}
    )
    if integrator.opts.adaptive
        zx_band_check(cache.β, cache.γ, cache.evalcache.αm, cache.evalcache.αf)
    end
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(
        integrator, cache::Union{NewmarkBetaCache, GeneralizedAlphaCache},
        repeat_step = false
    )
    (; t, dt, u, f, p) = integrator
    (; β, γ, thread, nlcache, evalcache, atmp) = cache

    # Evaluate predictor
    vₙ, uₙ = integrator.uprev.x
    if integrator.opts.adaptive
        # After a rejected attempt fsalfirst holds the trial state's derivative
        # (fsalfirst and fsallast share a buffer), so the left-endpoint
        # acceleration is always refreshed from uprev.
        f(integrator.fsalfirst, integrator.uprev, p, t)
        integrator.stats.nf += 1
    else
        # Only reachable with adaptive stepping off; the old guard
        # `derivative_discontinuity || !adaptive` was constant-true here and
        # reduces to `else`. Fixed-step keeps the historical refresh at the
        # trial state, pinned bit-for-bit by the fixed-step anchor tests.
        f(integrator.fsalfirst, integrator.u, p, t + dt)
        integrator.stats.nf += 1
    end
    aₙ = integrator.fsalfirst.x[1]

    # The AD buffers inside `evalcache` have to survive across steps, and not every inner
    # cache type exposes the parameter object it was handed, so the discretization state
    # is owned here and refreshed in place rather than rebuilt out of `nlcache`.
    evalcache.f = f
    evalcache.t = t
    evalcache.p = p
    evalcache.dt = dt
    evalcache.aₙ = aₙ
    evalcache.vₙ = vₙ
    evalcache.uₙ = uₙ

    SciMLBase.reinit!(nlcache, aₙ, p = evalcache)
    # A polyalgorithm cache delegates to its subcaches and a no-init cache delegates to a
    # one-shot `solve`, so neither holds the accepted iterate and both leave their own
    # retcode at `Default`; the outcome only exists on the solution `solve!` returns.
    nlsol = solve!(nlcache)
    if nlsol.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlsol.u

    @.. thread = thread u.x[1] = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)
    @.. thread = thread u.x[2] = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)

    if integrator.opts.adaptive
        κ = zx_kappa(evalcache.αm, evalcache.αf)
        Cp = dt^2 * zx_position_coefficient(β, κ)
        @.. thread = thread atmp.x[1] = Cp * (aₙ₊₁ - aₙ)
        calculate_residuals!(
            atmp.x[2], atmp.x[1], uₙ, u.x[2],
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        est = integrator.opts.internalnorm(atmp.x[2], t)
        # Velocity channel, clamped at the band edge so it stays live on the
        # family γ = 1/(2κ) where the measured coefficient vanishes; see the
        # derivation at the top of this file.
        Cv = dt * max(abs(zx_velocity_coefficient(γ, κ)), ZX_VELOCITY_BAND)
        @.. thread = thread atmp.x[1] = Cv * (aₙ₊₁ - aₙ)
        calculate_residuals!(
            atmp.x[2], atmp.x[1], vₙ, u.x[1],
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        est = max(est, integrator.opts.internalnorm(atmp.x[2], t))
        # A vanishing acceleration difference carries no step-size information, and
        # EEst = 0 makes the controller grow dt by qmax every step; a marginal
        # accept holds dt inside the controller's qsteady window instead.
        iszero(est) && (est = one(est))
        OrdinaryDiffEqCore.set_EEst!(integrator, est)
    end

    return
end

function initialize!(integrator, cache::NewmarkBetaConstantCache)
    if integrator.opts.adaptive
        zx_band_check(cache.β, cache.γ, zero(cache.β), zero(cache.β))
    end
    cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(
        integrator, cache::NewmarkBetaConstantCache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    (; β, γ, thread, nlsolver) = cache

    # Evaluate predictor
    if integrator.opts.adaptive
        integrator.fsalfirst .= f(integrator.uprev, p, t)
        integrator.stats.nf += 1
    else
        # Only reachable with adaptive stepping off; the old guard
        # `derivative_discontinuity || !adaptive` was constant-true here and
        # reduces to `else`. Fixed-step keeps the historical refresh at the
        # trial state, pinned bit-for-bit by the fixed-step anchor tests.
        integrator.fsalfirst .= f(integrator.u, p, t + dt)
        integrator.stats.nf += 1
    end
    aₙ = integrator.fsalfirst.x[1] # = fsallast
    vₙ, uₙ = integrator.uprev.x

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        zero(β), zero(β),
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )
    prob = NonlinearProblem{false}(discretized_residual, aₙ, evalcache)
    nlsol = solve(prob, nlsolver)
    if nlsol.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlsol.u

    # `update_uprev!` rebinds uprev to u for out-of-place problems instead of
    # copying, so the step must build fresh arrays and rebind u; mutating through
    # the old binding would drag uprev along and collapse the interpolation
    # interval to two identical endpoints.
    vₙ₊₁ = similar(vₙ)
    uₙ₊₁ = similar(uₙ)
    @.. thread = thread uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. thread = thread vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)
    integrator.u = ArrayPartition((vₙ₊₁, uₙ₊₁))

    if integrator.opts.adaptive
        zx_constant_cache_estimate!(
            integrator, β, γ, zero(β), zero(β), dt, t, aₙ, aₙ₊₁, vₙ, uₙ, vₙ₊₁, uₙ₊₁
        )
    end

    return
end

# Shared error estimate for the out-of-place steps; see the coefficient
# derivation at the top of this file.
function zx_constant_cache_estimate!(
        integrator, β, γ, αm, αf, dt, t, aₙ, aₙ₊₁, vₙ, uₙ, vₙ₊₁, uₙ₊₁
    )
    κ = zx_kappa(αm, αf)
    Cp = dt^2 * zx_position_coefficient(β, κ)
    tmp = similar(uₙ₊₁)
    @.. tmp = Cp * (aₙ₊₁ - aₙ)
    atmp = calculate_residuals(
        tmp, uₙ, uₙ₊₁,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )
    est = integrator.opts.internalnorm(atmp, t)
    # Velocity channel, clamped at the band edge so it stays live on the
    # family γ = 1/(2κ) where the measured coefficient vanishes; see the
    # derivation at the top of this file.
    Cv = dt * max(abs(zx_velocity_coefficient(γ, κ)), ZX_VELOCITY_BAND)
    @.. tmp = Cv * (aₙ₊₁ - aₙ)
    atmp = calculate_residuals(
        tmp, vₙ, vₙ₊₁,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )
    est = max(est, integrator.opts.internalnorm(atmp, t))
    # A vanishing acceleration difference carries no step-size information, and
    # EEst = 0 makes the controller grow dt by qmax every step; a marginal
    # accept holds dt inside the controller's qsteady window instead.
    iszero(est) && (est = one(est))
    OrdinaryDiffEqCore.set_EEst!(integrator, est)
    return nothing
end

function initialize!(integrator, cache::GeneralizedAlphaConstantCache)
    if integrator.opts.adaptive
        zx_band_check(cache.β, cache.γ, cache.αm, cache.αf)
    end
    cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(
        integrator, cache::GeneralizedAlphaConstantCache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    (; αm, αf, β, γ, thread, nlsolver) = cache

    if integrator.opts.adaptive
        integrator.fsalfirst .= f(integrator.uprev, p, t)
        integrator.stats.nf += 1
    else
        # Only reachable with adaptive stepping off; the old guard
        # `derivative_discontinuity || !adaptive` was constant-true here and
        # reduces to `else`. Fixed-step keeps the historical refresh at the
        # trial state, pinned bit-for-bit by the fixed-step anchor tests.
        integrator.fsalfirst .= f(integrator.u, p, t + dt)
        integrator.stats.nf += 1
    end
    aₙ = integrator.fsalfirst.x[1]
    vₙ, uₙ = integrator.uprev.x

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        αm, αf,
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )
    prob = NonlinearProblem{false}(discretized_residual, aₙ, evalcache)
    nlsol = solve(prob, nlsolver)
    if nlsol.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlsol.u

    # `update_uprev!` rebinds uprev to u for out-of-place problems instead of
    # copying, so the step must build fresh arrays and rebind u; mutating through
    # the old binding would drag uprev along and collapse the interpolation
    # interval to two identical endpoints.
    vₙ₊₁ = similar(vₙ)
    uₙ₊₁ = similar(uₙ)
    @.. thread = thread uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. thread = thread vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)
    integrator.u = ArrayPartition((vₙ₊₁, uₙ₊₁))

    if integrator.opts.adaptive
        zx_constant_cache_estimate!(
            integrator, β, γ, αm, αf, dt, t, aₙ, aₙ₊₁, vₙ, uₙ, vₙ₊₁, uₙ₊₁
        )
    end

    return
end

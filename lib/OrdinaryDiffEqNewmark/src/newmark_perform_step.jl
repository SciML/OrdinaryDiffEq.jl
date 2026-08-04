# Newmark-β is generalized-α with αm = αf = 0, and both schemes reach the inner solver
# through the same discretization cache, so the in-place step is shared.
function initialize!(
        integrator, cache::Union{NewmarkBetaCache, GeneralizedAlphaCache}
    )
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
    if integrator.derivative_discontinuity || !integrator.opts.adaptive
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
        f(integrator.fsallast, u, p, t + dt)
        integrator.stats.nf += 1

        if integrator.success_iter == 0
            OrdinaryDiffEqCore.set_EEst!(integrator, one(OrdinaryDiffEqCore.get_EEst(integrator)))
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            @.. thread = thread atmp = (integrator.fsallast.x[1] - aₙ₊₁)
            OrdinaryDiffEqCore.set_EEst!(
                integrator,
                dt * dt * (β - 1 // 6) *
                    integrator.opts.internalnorm(atmp, t)
            )
        end
    end

    return
end

function initialize!(integrator, cache::NewmarkBetaConstantCache)
    cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(
        integrator, cache::NewmarkBetaConstantCache, repeat_step = false
    )
    (; t, u, dt, f, p) = integrator
    (; β, γ, thread, nlsolver, atmp) = cache

    # Evaluate predictor
    if integrator.derivative_discontinuity || !integrator.opts.adaptive
        integrator.fsalfirst .= f(u, p, t + dt)
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

    # The velocity component in uprev and u is shadowed, so the order of these two operation below matter.
    @.. thread = thread u.x[2] = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. thread = thread u.x[1] = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    if integrator.opts.adaptive
        integrator.fsallast .= f(u, p, t + dt)
        integrator.stats.nf += 1

        if integrator.success_iter == 0
            OrdinaryDiffEqCore.set_EEst!(integrator, one(OrdinaryDiffEqCore.get_EEst(integrator)))
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            @.. thread = thread atmp = (integrator.fsallast.x[1] - aₙ₊₁)
            OrdinaryDiffEqCore.set_EEst!(
                integrator,
                dt * dt * (β - 1 // 6) *
                    integrator.opts.internalnorm(atmp, t)
            )
        end
    end

    return
end

function initialize!(integrator, cache::GeneralizedAlphaConstantCache)
    cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(
        integrator, cache::GeneralizedAlphaConstantCache, repeat_step = false
    )
    (; t, u, dt, f, p) = integrator
    (; αm, αf, β, γ, thread, nlsolver, atmp) = cache

    if integrator.derivative_discontinuity || !integrator.opts.adaptive
        integrator.fsalfirst .= f(u, p, t + dt)
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

    # velocity component in uprev and u is shadowed, so order matters here
    @.. thread = thread u.x[2] = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. thread = thread u.x[1] = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    if integrator.opts.adaptive
        integrator.fsallast .= f(u, p, t + dt)
        integrator.stats.nf += 1

        if integrator.success_iter == 0
            OrdinaryDiffEqCore.set_EEst!(integrator, one(OrdinaryDiffEqCore.get_EEst(integrator)))
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            @.. thread = thread atmp = (integrator.fsallast.x[1] - aₙ₊₁)
            OrdinaryDiffEqCore.set_EEst!(
                integrator,
                dt * dt * (β - 1 // 6) *
                    integrator.opts.internalnorm(atmp, t)
            )
        end
    end

    return
end

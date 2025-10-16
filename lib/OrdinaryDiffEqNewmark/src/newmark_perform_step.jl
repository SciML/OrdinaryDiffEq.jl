function initialize!(integrator, cache::NewmarkBetaCache)
    duprev, uprev = integrator.uprev.x
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    # integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    return
end

@muladd function perform_step!(integrator, cache::NewmarkBetaCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    @unpack β, γ, nlcache = cache

    M = f.mass_matrix

    # Evaluate predictor
    if integrator.u_modified
        f(integrator.fsalfirst.x[1], u, p, t + dt)
    end
    aₙ = integrator.fsalfirst.x[1]
    vₙ, uₙ = integrator.uprev.x

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ
    )
    SciMLBase.reinit!(nlcache, aₙ, p = evalcache)
    solve!(nlcache)
    if nlcache.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlcache.u

    @.. integrator.u = ArrayPartition(
        vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁),
        uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    )

    f(integrator.fsallast, u, p, t + dt)
    integrator.stats.nf += 1

    #
    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            δaₙ₊₁ = (integrator.fsallast.x[1] - aₙ₊₁)
            integrator.EEst = dt * dt * (β - 1 // 6) *
                              integrator.opts.internalnorm(δaₙ₊₁, t)
        end
    end

    return
end

function initialize!(integrator, cache::NewmarkBetaConstantCache)
    duprev, uprev = integrator.uprev.x
    cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    # integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    return
end

@muladd function perform_step!(
        integrator, cache::NewmarkBetaConstantCache, repeat_step = false)
    @unpack t, u, dt, f, p = integrator
    @unpack β, γ, nlsolver = cache

    M = f.mass_matrix

    # Evaluate predictor
    if integrator.u_modified
        @.. integrator.fsalfirst.x[1] = f(u, p, t + dt)
    end
    aₙ = integrator.fsalfirst.x[1]
    vₙ, uₙ = integrator.uprev.x

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ
    )
    prob = NonlinearProblem{false}(newmark_discretized_residual, aₙ, evalcache)
    nlsol = solve(prob, nlsolver)
    if nlsol.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlsol.u

    @.. u = ArrayPartition(
        vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁),
        uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    )

    integrator.fsallast .= f(u, p, t + dt)
    integrator.stats.nf += 1

    #
    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            δaₙ₊₁ = (integrator.fsallast.x[1] - aₙ₊₁)
            integrator.EEst = dt * dt / 2 * (2 * β - 1 / 3) *
                              integrator.opts.internalnorm(δaₙ₊₁, t)
        end
    end

    return
end

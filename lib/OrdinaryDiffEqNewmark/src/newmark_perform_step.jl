function initialize!(integrator, cache::NewmarkBetaCache)
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.fsalfirst
    return
end

@muladd function perform_step!(integrator, cache::NewmarkBetaCache, repeat_step = false)
    (; t, dt, u, f, p) = integrator
    (; β, γ, thread, nlcache, atmp) = cache

    # Evaluate predictor
    vₙ, uₙ = integrator.uprev.x
    if integrator.u_modified || !integrator.opts.adaptive
        f(integrator.fsalfirst, integrator.u, p, t + dt)
        integrator.stats.nf += 1
    end
    aₙ = integrator.fsalfirst.x[1]

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
        nlcache.p.atmp, nlcache.p.vₙ₊₁,  nlcache.p.uₙ₊₁,
    )
    SciMLBase.reinit!(nlcache, aₙ, p = evalcache)
    solve!(nlcache)
    if nlcache.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlcache.u

    @.. thread=thread u.x[1] = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)
    @.. thread=thread u.x[2] = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)

    if integrator.opts.adaptive
        f(integrator.fsallast, u, p, t + dt)
        integrator.stats.nf += 1

        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            @.. thread=thread atmp = (integrator.fsallast - aₙ₊₁)
            integrator.EEst = dt * dt * (β - 1 // 6) *
                              integrator.opts.internalnorm(atmp, t)
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
        integrator, cache::NewmarkBetaConstantCache, repeat_step = false)
    (; t, u, dt, f, p) = integrator
    (; β, γ, thread, nlsolver, atmp) = cache

    # Evaluate predictor
    if integrator.u_modified || !integrator.opts.adaptive
        integrator.fsalfirst .= f(u, p, t + dt)
        integrator.stats.nf += 1
    end
    aₙ = integrator.fsalfirst.x[1] # = fsallast
    vₙ, uₙ = integrator.uprev.x

    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )
    prob = NonlinearProblem{false}(newmark_discretized_residual, aₙ, evalcache)
    nlsol = solve(prob, nlsolver)
    if nlsol.retcode != ReturnCode.Success
        integrator.force_stepfail = true
        return
    end
    aₙ₊₁ = nlsol.u

    # The velocity component in uprev and u is shadowed, so the order of these two operation below matter.
    @.. thread=thread u.x[2] = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. thread=thread u.x[1] = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    # @info "A", integrator.opts.internalnorm(aₙ₊₁,t), integrator.opts.internalnorm(vₙ,t), integrator.opts.internalnorm(aₙ,t)
    # u .= ArrayPartition(
    #     vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁),
    #     uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁),
    # )
    # @info integrator.opts.internalnorm(aₙ₊₁,t), integrator.opts.internalnorm(vₙ,t), integrator.opts.internalnorm(aₙ,t)

    # # u.x[1] .= vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)
    # # u.x[2] .= uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)

    #
    if integrator.opts.adaptive
        integrator.fsallast .= f(u, p, t + dt)
        integrator.stats.nf += 1

        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            @.. thread=thread atmp = (integrator.fsallast.x[1] - aₙ₊₁)
            integrator.EEst = dt * dt * (β - 1 // 6) *
                              integrator.opts.internalnorm(atmp, t)
        end
    end

    return
end

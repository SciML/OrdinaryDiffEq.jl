function initialize!(integrator, cache::NewmarkBetaCache)
    duprev, uprev = integrator.uprev.x
    if isinplace(integrator.f)
        integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    else
        cache.fsalfirst .= integrator.f(integrator.uprev, integrator.p, integrator.t)
    end
    integrator.stats.nf += 1
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast  = cache.fsalfirst
    # integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    return
end

@muladd function perform_step!(integrator, cache::NewmarkBetaCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    @unpack upred, β, γ, nlsolver = cache

    M = f.mass_matrix

    # Evaluate predictor
    aₙ     = integrator.fsalfirst.x[1]
    vₙ, uₙ = integrator.uprev.x

    # Manually unrolled to see what needs to go where
    aₙ₊₁ = copy(aₙ) # acceleration term
    atmp = copy(aₙ)

    nlf = isinplace(f) ? newmark_discretized_residual! : newmark_discretized_residual
    nlprob = NonlinearProblem{isinplace(f)}(nlf, aₙ, NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
        atmp, copy(uₙ), copy(vₙ),
    ))
    nlsol = solve(nlprob, nlsolver)
    aₙ₊₁ = nlsol.u

    u = ArrayPartition(
        vₙ + dt * ((1-γ)*aₙ + γ*aₙ₊₁),
        uₙ + dt * vₙ + dt^2/2 * ((1-2β)*aₙ + 2β*aₙ₊₁),
    )

    if isinplace(f)
        f(integrator.fsallast, u, p, t + dt)
    else
        integrator.fsallast .= f(u, p, t + dt)
    end
    integrator.stats.nf += 1
    integrator.u = u

    #
    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            # Zienkiewicz and Xie (1991) Eq. 21
            δaₙ₊₁ = (integrator.fsallast.x[1] - aₙ₊₁)
            integrator.EEst = dt*dt/2 * (2*β - 1/3) * integrator.opts.internalnorm(δaₙ₊₁, t)
        end
    end

    return
end

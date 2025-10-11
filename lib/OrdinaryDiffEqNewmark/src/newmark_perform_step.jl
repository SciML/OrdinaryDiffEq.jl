function initialize!(integrator, cache::NewmarkBetaCache)
    duprev, uprev = integrator.uprev.x
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
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
    # J = zeros(length(aₙ), length(aₙ))
    # for i in 1:10 # = max iter - Newton loop for eq [1] above
    #     uₙ₊₁ = uₙ + dt * vₙ + dt^2/2 * ((1-2β)*aₙ + 2β*aₙ₊₁)
    #     vₙ₊₁ = vₙ + dt * ((1-γ)*aₙ + γ*aₙ₊₁)
    #     # Compute residual
    #     f.f1(atmp, vₙ₊₁, uₙ₊₁, p, t)
    #     integrator.stats.nf += 1
    #     residual = M*(aₙ₊₁ - atmp)
    #     # Compute jacobian
    #     f.jac(J, vₙ₊₁, uₙ₊₁, (γ*dt, β*dt*dt), p, t)
    #     # Solve for increment
    #     Δaₙ₊₁ = (M-J) \ residual
    #     aₙ₊₁ .-= Δaₙ₊₁ # Looks like I messed up the signs somewhere :')
    #     increment_norm = integrator.opts.internalnorm(Δaₙ₊₁, t)
    #     increment_norm < 1e-4 && break
    #     i == 10 && error("Newton diverged. ||Δaₙ₊₁||=$increment_norm")
    # end

    nlf = isinplace(f) ? newmark_discretized_residual! : newmark_discretized_residual
    nlprob = NonlinearProblem{isinplace(f)}(nlf, aₙ, NewmarkDiscretizationCache(;
        f, dt, t, p,
        vₙ, uₙ, aₙ,
        uₙ₊₁ = copy(uₙ),
        vₙ₊₁ = copy(vₙ),
        β, γ,
        atmp,
    ))
    nlsol = solve(nlprob, nlsolver)
    aₙ₊₁ = nlsol.u

    u = ArrayPartition(
        vₙ + dt * ((1-γ)*aₙ + γ*aₙ₊₁),
        uₙ + dt * vₙ + dt^2/2 * ((1-2β)*aₙ + 2β*aₙ₊₁),
    )

    f(integrator.fsallast, u, p, t + dt)
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

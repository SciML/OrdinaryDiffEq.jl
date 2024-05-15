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

    # Given Mu'' = f(u,u',t) we need to solve the non-linear problem
    #   Mu(tₙ₊₁)'' - f(ũ(tₙ₊₁) + u(tₙ₊₁)'' β Δtₙ²,ũ(tₙ₊₁)' + u(tₙ₊₁)'' γ Δtₙ,t) = 0
    # for u(tₙ₊₁)''.
    # The predictors ũ(tₙ₊₁) are computed as
    #   ũ(tₙ₊₁)  = u(tₙ) + u(tₙ)' Δtₙ + u(tₙ)'' (0.5 - β) Δtₙ²
    #   ũ(tₙ₊₁)' = u(tₙ)' + u(tₙ)'' (1.0 - γ) Δtₙ
    # such that we can compute the solution with the correctors
    #   u(tₙ₊₁)  = ũ(tₙ₊₁)  + u(tₙ₊₁)'' β Δtₙ²
    #   u(tₙ₊₁)' = ũ(tₙ₊₁)' + u(tₙ₊₁)'' γ Δtₙ

    mass_matrix = f.mass_matrix

    # Evaluate predictor
    dduprev = integrator.fsalfirst.x[1]
    duprev, uprev = integrator.uprev.x
    upred_full = ArrayPartition(
        duprev + dt*(1.0 - γ)*dduprev,
        uprev  + dt*dt*(0.5 - β)*dduprev + dt*duprev
    )

    # _tmp = mass_matrix * @.. broadcast=false (α₁ * uprev+α₂ * uprev2)
    # nlsolver.tmp = @.. broadcast=false _tmp/(dt * β₀)

    # nlsolve!(...) solves for
    #   dt⋅f(innertmp + γ̂⋅z, p, t + c⋅dt) + outertmp = z
    # So we rewrite the problem
    #     u(tₙ₊₁)'' - f₁(ũ(tₙ₊₁) + u(tₙ₊₁)'' β Δtₙ², ũ(tₙ₊₁)' + u(tₙ₊₁)'' γ Δtₙ,t) = 0
    #   z = Δtₙ u(tₙ₊₁)'':
    #     z         - Δtₙ f₁(ũ(tₙ₊₁) +         z β Δtₙ, ũ(tₙ₊₁)' +         z γ,t) = 0
    #                 Δtₙ f₁(ũ(tₙ₊₁) +         z β Δtₙ, ũ(tₙ₊₁)' +         z γ,t) = z
    #   γ̂ = [γ, β Δtₙ]:
    #                 Δtₙ f₁(ũ(tₙ₊₁) +         z γ̂₂    , ũ(tₙ₊₁)' +         z γ̂₁   ,t) = z
    #   innertmp = [ũ(tₙ₊₁)', ũ(tₙ₊₁)]:
    #                 Δtₙ f₁(innertmp₂ +         z β Δtₙ², innertmp₁ +         z γ Δtₙ,t) = z
    # Note: innertmp = nlsolve.tmp
    nlsolver.γ = ArrayPartitionNLSolveHelper(γ, β * dt) # = γ̂
    # nlsolver.γ = ArrayPartitionNLSolveHelper(0.0, β * dt)
    nlsolver.tmp .= upred_full # TODO check f tmp is potentially modified and if not elimiate the allocation of upred_full

    # Use the linear extrapolation Δtₙ u(tₙ)'' as initial guess for the nonlinear solve
    nlsolver.z = dt*dduprev
    ddu = nlsolve!(nlsolver, integrator, cache, repeat_step) / dt
    nlsolvefail(nlsolver) && return

    # Apply corrector
    u = ArrayPartition(
        upred_full.x[1] + ddu*γ*dt,
        upred_full.x[2] + ddu*β*dt*dt,
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
            δddu = (integrator.fsallast.x[1] - ddu)
            integrator.EEst = dt*dt/2 * (2*β - 1/3) * integrator.opts.internalnorm(δddu, t)
        end
    end

    return
end

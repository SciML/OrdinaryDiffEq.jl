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

    # This is derived from the idea stated in Nonlinear Finite Elements by Peter Wriggers, Ch 6.1.2 .
    #
    # Let us introduce the notation v = u' and a = u'' = v' such that we write the ODE problem as Ma = f(u,v,t).
    # For the time discretization we assume that:
    #   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 aₙ₊ₐ₁
    #   vₙ₊₁ = vₙ + Δtₙ aₙ₊ₐ₂
    # with a₁ = 1-2β and a₂ = 1-γ, such that
    #   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁]
    #   vₙ₊₁ = vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁]
    # 
    # This allows us to reduce the implicit discretization to have only aₙ₊₁ as the unknown:
    #   Maₙ₊₁ = f(uₙ₊₁(aₙ₊₁), vₙ₊₁(aₙ₊₁), tₙ₊₁) 
    #         = f(uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁], vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁], tₙ₊₁)
    # Such that we have to solve the nonlinear problem
    #   Maₙ₊₁ - f(uₙ₊₁(aₙ₊₁), vₙ₊₁(aₙ₊₁), tₙ₊₁)  = 0
    # for aₙ₊₁'' in each time step.

    # For the Newton method the linearization becomes
    #   M - (dₐuₙ₊₁ ∂fᵤ + dₐvₙ₊₁ ∂fᵥ) = 0
    #   M - (Δtₙ²β  ∂fᵤ +  Δtₙγ  ∂fᵥ) = 0

    M = f.mass_matrix

    # Evaluate predictor
    aₙ     = integrator.fsalfirst.x[1]
    uₙ, vₙ = integrator.uprev.x

    # _tmp = mass_matrix * @.. broadcast=false (α₁ * uprev+α₂ * uprev2)
    # nlsolver.tmp = @.. broadcast=false _tmp/(dt * β₀)

    # Note, we switch to notation closer to the SciML implemenation now. Needs to be double checked, also to be consistent with the formulation above
    # nlsolve!(...) solves for
    #   dt⋅f(innertmp + γ̂⋅z, p, t + c⋅dt) + outertmp = z
    # So we rewrite the problem
    #     u(tₙ₊₁)'' - f₁(ũ(tₙ₊₁) + u(tₙ₊₁)'' 2β Δtₙ², ũ(tₙ₊₁)' + u(tₙ₊₁)'' γ Δtₙ,t) = 0
    #   z = Δtₙ u(tₙ₊₁)'':
    #     z         - Δtₙ f₁(ũ(tₙ₊₁) +         z 2β Δtₙ, ũ(tₙ₊₁)' +         z γ,t) = 0
    #                 Δtₙ f₁(ũ(tₙ₊₁) +         z 2β Δtₙ, ũ(tₙ₊₁)' +         z γ,t) = z
    #   γ̂ = [γ, 2β Δtₙ]:
    #                 Δtₙ f₁(ũ(tₙ₊₁) +         z γ̂₂    , ũ(tₙ₊₁)' +         z γ̂₁   ,t) = z
    #   innertmp = [ũ(tₙ₊₁)', ũ(tₙ₊₁)]:
    #                 Δtₙ f₁(innertmp₂ +         z 2β Δtₙ², innertmp₁ +         z γ Δtₙ,t) = z
    # Note: innertmp = nlsolve.tmp
    # nlsolver.γ = ???
    # nlsolver.tmp .= vₙ # TODO check f tmp is potentially modified and if not elimiate the allocation of upred_full
    # nlsolver.z .= aₙ
    # ddu = nlsolve!(nlsolver, integrator, cache, repeat_step) / dt
    # nlsolvefail(nlsolver) && return

    # Manually unrolled to see what needs to go where
    aₙ₊₁ = copy(aₙ) # acceleration term
    atmp = copy(aₙ)
    for i in 1:10 # = max iter - Newton loop for eq [1] above
        uₙ₊₁ = uₙ + Δtₙ * vₙ + Δtₙ^2/2 * ((1-2β)*aₙ + 2β*aₙ₊₁)
        vₙ₊₁ = vₙ + Δtₙ * ((1-γ)aₙ + γaₙ₊₁)
        # Compute residual
        f.f1(atmp, uₙ₊₁, vₙ₊₁, p, t)
        residual = M*(aₙ₊₁ - atmp)
        # Compute jacobian
        f.jac(J, upred_full, (β*dt*dt, γ*dt), t)
        # Solve for increment
        Δddu = (M-J) \ residual
        ddu .+= Δddu
        norm(Δddu) < 1e-4 && break
        i == 10 && error("Newton diverged. Δddu=$(norm(Δddu))")
    end
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

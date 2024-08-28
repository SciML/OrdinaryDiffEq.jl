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
    # Let us introduce the notation v = u' and a = u'' = v' such that we write the ODE problem as Ma = f(v,u,t).
    # For the time discretization we assume that:
    #   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 aₙ₊ₐ₁
    #   vₙ₊₁ = vₙ + Δtₙ aₙ₊ₐ₂
    # with a₁ = 1-2β and a₂ = 1-γ, such that
    #   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁]
    #   vₙ₊₁ = vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁]
    # 
    # This allows us to reduce the implicit discretization to have only aₙ₊₁ as the unknown:
    #   Maₙ₊₁ = f(vₙ₊₁(aₙ₊₁), uₙ₊₁(aₙ₊₁), tₙ₊₁) 
    #         = f(vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁], uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁], tₙ₊₁)
    # Such that we have to solve the nonlinear problem
    #   Maₙ₊₁ - f(vₙ₊₁(aₙ₊₁), uₙ₊₁(aₙ₊₁), tₙ₊₁)  = 0
    # for aₙ₊₁'' in each time step.

    # For the Newton method the linearization becomes
    #   M - (dₐuₙ₊₁ ∂fᵤ + dₐvₙ₊₁ ∂fᵥ) = 0
    #   M - (Δtₙ²β  ∂fᵤ +  Δtₙγ  ∂fᵥ) = 0

    M = f.mass_matrix

    # Evaluate predictor
    aₙ     = integrator.fsalfirst.x[1]
    vₙ, uₙ = integrator.uprev.x

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
    # aₙ₊₁ = nlsolve!(nlsolver, integrator, cache, repeat_step) / dt
    # nlsolvefail(nlsolver) && return

    # Manually unrolled to see what needs to go where
    aₙ₊₁ = copy(aₙ) # acceleration term
    atmp = copy(aₙ)
    J = zeros(length(aₙ), length(aₙ))
    for i in 1:10 # = max iter - Newton loop for eq [1] above
        uₙ₊₁ = uₙ + dt * vₙ + dt^2/2 * ((1-2β)*aₙ + 2β*aₙ₊₁)
        vₙ₊₁ = vₙ + dt * ((1-γ)*aₙ + γ*aₙ₊₁)
        # Compute residual
        f.f1(atmp, vₙ₊₁, uₙ₊₁, p, t)
        integrator.stats.nf += 1
        residual = M*(aₙ₊₁ - atmp)
        # Compute jacobian
        f.jac(J, vₙ₊₁, uₙ₊₁, (γ*dt, β*dt*dt), p, t)
        # Solve for increment
        Δaₙ₊₁ = (M-J) \ residual
        aₙ₊₁ .-= Δaₙ₊₁ # Looks like I messed up the signs somewhere :')
        increment_norm = integrator.opts.internalnorm(Δaₙ₊₁, t)
        increment_norm < 1e-4 && break
        i == 10 && error("Newton diverged. ||Δaₙ₊₁||=$increment_norm")
    end

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

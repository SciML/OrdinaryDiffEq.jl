@muladd function perform_step!(integrator, cache::SKenCarpConstantCache)
    (; t, dt, uprev, u, p, f) = integrator
    g = integrator.f.g
    (; nlsolver) = cache
    (;
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32,
    ) = cache.tab
    (;
        ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4,
        ebtilde1, ebtilde2, ebtilde3, ebtilde4,
    ) = cache.tab
    (; nb021, nb043) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    sqrt3 = sqrt(3one(eltype(integrator.W.dW)))
    chi2 = (integrator.W.dW + integrator.W.dZ / sqrt3) / 2 #I_(1,0)/h

    if integrator.f isa SplitSDEFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    # calculate W
    repeat_step = false

    z₁ = dt * f(uprev, p, t)
    nlsolver.c = 2γ

    g1 = g(uprev, p, t)
    tmp = uprev + γ * z₁ + nb021 * chi2 .* g1
    if integrator.f isa SplitSDEFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * f2(uprev, p, t)
        tmp += ea21 * k1
    end
    nlsolver.tmp = tmp

    if alg.extrapolant == :min_correct
        z₂ = zero(z₁)
    elseif alg.extrapolant == :trivial
        z₂ = z₁
    end
    nlsolver.z = z₂

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.c = c3
    if integrator.f isa SplitSDEFunction
        u = tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γ * dt)
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.tmp = tmp

    if alg.extrapolant == :min_correct
        z₃ = zero(z₂)
    elseif alg.extrapolant == :trivial
        z₃ = z₂
    end
    nlsolver.z = z₃

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.c = one(nlsolver.c)

    # Note: Can use g1 since c13 = 0 => g3 == g1

    if integrator.f isa SplitSDEFunction
        u = tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3 +
            nb043 * chi2 .* g1
    else
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + nb043 * chi2 .* g1
    end
    nlsolver.tmp = tmp

    if alg.extrapolant == :min_correct
        z₄ = zero(z₂)
    elseif alg.extrapolant == :trivial
        z₄ = z₂
    end
    nlsolver.z = z₄

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = tmp + γ * z₄
    g4 = g(uprev, p, t + dt)

    E₂ = chi2 .* (g1 - g4)

    if integrator.f isa SplitSDEFunction
        k4 = dt * f2(u, p, t + dt)
        u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 + eb2 * k2 + eb3 * k3 + eb4 * k4 +
            integrator.W.dW .* g4 + E₂
    else
        u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + integrator.W.dW .* g4 + E₂
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if alg.ode_error_est
            if integrator.f isa SplitSDEFunction
                tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + ebtilde1 * k1 +
                    ebtilde2 * k2 + ebtilde3 * k3 + ebtilde4 * k4 + chi2 * (g1 - g4)
            else
                tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + chi2 * (g1 - g4)
            end
            if alg.smooth_est # From Shampine
                E₁ = DiffEqBase._reshape(
                    get_W(nlsolver) \ DiffEqBase._vec(tmp),
                    axes(tmp)
                )
            else
                E₁ = tmp
            end
        else
            E₁ = z₁ + z₂ + z₃ + z₄
        end

        resids = calculate_residuals(
            E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(resids, t))
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SKenCarpCache)
    (; t, dt, uprev, u, p, f) = integrator
    g = integrator.f.g
    (; z₁, z₂, z₃, z₄, k1, k2, k3, k4, atmp) = cache
    (; g1, g4, chi2, nlsolver) = cache
    (; z, tmp) = nlsolver
    (; k, dz) = nlsolver.cache # alias to reduce memory
    (;
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32,
    ) = cache.tab
    (; ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4) = cache.tab
    (; ebtilde1, ebtilde2, ebtilde3, ebtilde4) = cache.tab
    (; nb021, nb043) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    repeat_step = false

    # Some aliases

    E₁ = g4
    E₂ = dz

    if integrator.f isa SplitSDEFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    sqrt3 = sqrt(3one(eltype(integrator.W.dW)))
    if integrator.W.dW isa Union{SArray, Number}
        chi2 = (integrator.W.dW + integrator.W.dZ / sqrt3) / 2 #I_(1,0)/h
    else
        @.. chi2 = (integrator.W.dW + integrator.W.dZ / sqrt3) / 2 #I_(1,0)/h
    end

    # precalculations

    γdt = γ * dt

    if !repeat_step && !integrator.last_stepfail
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    end

    ##### Step 2

    # TODO: Add a cache so this isn't overwritten near the end, so it can not repeat on fail
    g(g1, uprev, p, t)

    if is_diagonal_noise(integrator.sol.prob)
        @.. z₄ = chi2 * g1 # use z₄ as storage for the g1*chi2
    else
        mul!(z₄, g1, chi2) # use z₄ as storage for the g1*chi2
    end

    @.. tmp = uprev + γ * z₁ + nb021 * z₄

    if alg.extrapolant == :min_correct
        @.. z₂ = zero(eltype(dz))
    elseif alg.extrapolant == :trivial
        @.. z₂ = z₁
    end

    if integrator.f isa SplitSDEFunction
        # This assumes the implicit part is cheaper than the explicit part
        if !repeat_step && !integrator.last_stepfail
            f2(k1, integrator.uprev, integrator.p, integrator.t)
            k1 .*= dt
        end
        @.. tmp += ea21 * k1
    end

    nlsolver.z = z₂
    nlsolver.c = 2γ

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitSDEFunction
        @.. u = tmp + γ * z₂
        f2(k2, u, p, t + 2γ * dt)
        k2 .*= dt
        @.. tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        @.. tmp = uprev + a31 * z₁ + a32 * z₂
    end

    if alg.extrapolant == :min_correct
        @.. z₃ = zero(eltype(dz))
    elseif alg.extrapolant == :trivial
        @.. z₃ = z₂
    end
    nlsolver.z = z₃
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitSDEFunction
        @.. u = tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        # z₄ is storage for the g1*chi2 from earlier
        @.. tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3 + nb043 * z₄
    else
        (; α41, α42) = cache.tab
        # z₄ is storage for the g1*chi2
        @.. tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + nb043 * z₄
    end

    if alg.extrapolant == :min_correct
        @.. z₄ = zero(eltype(dz))
    elseif alg.extrapolant == :trivial
        @.. z₄ = z₂
    end

    nlsolver.z = z₄
    nlsolver.c = one(nlsolver.c)
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    g(g4, u, p, t + dt)

    if integrator.f isa SplitSDEFunction
        @.. u = tmp + γ * z₄
        f2(k4, u, p, t + dt)
        k4 .*= dt
        if is_diagonal_noise(integrator.sol.prob)
            @.. E₂ = chi2 * (g1 - g4)
            @.. u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 + eb2 * k2 + eb3 * k3 +
                eb4 * k4 + integrator.W.dW * g4 + E₂
        else
            g1 .-= g4
            mul!(E₂, g1, chi2)
            mul!(tmp, g4, integrator.W.dW)
            @.. u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 + eb2 * k2 + eb3 * k3 +
                eb4 * k4 + tmp + E₂
        end
    else
        if is_diagonal_noise(integrator.sol.prob)
            @.. E₂ = chi2 * (g1 - g4)
            @.. u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + integrator.W.dW * g4 + E₂
        else
            g1 .-= g4
            mul!(E₂, g1, chi2)
            mul!(tmp, g4, integrator.W.dW)
            @.. u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + tmp + E₂
        end
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if alg.ode_error_est
            if integrator.f isa SplitSDEFunction
                @.. g1 = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + ebtilde1 * k1 +
                    ebtilde2 * k2 + ebtilde3 * k3 + ebtilde4 * k4
            else
                @.. g1 = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄
            end
            if alg.smooth_est # From Shampine
                linres = dolinsolve(
                    integrator, nlsolver.cache.linsolve;
                    b = DiffEqBase._vec(g1),
                    linu = DiffEqBase._vec(E₁)
                )
            else
                E₁ .= dz
            end
        else
            @.. E₁ = z₁ + z₂ + z₃ + z₄
        end

        calculate_residuals!(
            tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end
end

@muladd function perform_step!(
        integrator, cache::CFNLIRK3ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, p) = integrator
    nlsolver = cache.nlsolver
    (; γ, a31, a32, a41, a42, a43, c2, c3, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4) = cache.tab
    alg = unwrap_alg(integrator, true)

    f2 = nothing
    k1 = nothing
    k2 = nothing
    k3 = nothing
    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    # calculate W
    markfirststage!(nlsolver)

    if integrator.f isa SplitFunction
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        z₁ = dt .* f(uprev, p, t)
    else
        # FSAL Step 1
        z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation for guess
    nlsolver.z = z₂ = z₁

    nlsolver.tmp = uprev

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt .* f2(uprev, p, t)
        nlsolver.tmp += ea21 * k1
    end

    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + c2 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        z₃ = z₂
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃
    nlsolver.tmp = tmp
    nlsolver.c = c3

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ = z₃
        u = nlsolver.tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3
    else
        z₄ = z₃
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄
    nlsolver.c = 1
    nlsolver.tmp = tmp

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₄
    if integrator.f isa SplitFunction
        k4 = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 + eb2 * k2 +
            eb3 * k3 + eb4 * k4
    end

    ################################### Finalize

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z₄ ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::CFNLIRK3Cache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; z₁, z₂, z₃, z₄, k1, k2, k3, k4, atmp, nlsolver) = cache
    (; tmp) = nlsolver
    (; γ, a31, a32, a41, a42, a43, c2, c3) = cache.tab
    (; ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4) = cache.tab

    alg = unwrap_alg(integrator, true)

    f2 = nothing
    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast = false z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation for guess
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast = false tmp = uprev

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast = false k1 = dt * integrator.fsalfirst - z₁
        @.. broadcast = false tmp += ea21 * k1
    end

    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast = false u = tmp + γ * z₂
        f2(k2, u, p, t + c2 * dt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast = false tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        @.. broadcast = false z₃ = z₂
        @.. broadcast = false tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₂
        @.. broadcast = false u = tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast = false tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
            ea42 * k2 + ea43 * k3
    else
        @.. broadcast = false z₄ = z₂
        @.. broadcast = false tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = 1
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast = false u = tmp + γ * z₄
    if integrator.f isa SplitFunction
        f2(k4, u, p, t + dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast = false u = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 +
            eb2 * k2 + eb3 * k3 + eb4 * k4
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast = false integrator.fsallast = z₄ / dt
    end
end

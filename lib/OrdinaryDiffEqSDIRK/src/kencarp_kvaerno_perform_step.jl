@muladd function perform_step!(integrator, cache::Kvaerno3ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31, α32 = cache.tab
    alg = unwrap_alg(integrator, true)

    # calculate W
    markfirststage!(nlsolver)

    # FSAL Step 1
    nlsolver.z = z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    nlsolver.z = z₂ = z₁

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess is from Hermite derivative on z₁ and z₂
    nlsolver.z = z₃ = α31 * z₁ + α32 * z₂

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = a31 * z₁ + a32 * z₂ + γ * z₃ # use yhat as prediction

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = 1
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₄

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₄ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Kvaerno3Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, z₂, z₃, z₄, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31, α32 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    # FSAL Step 1
    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    # Guess is from Hermite derivative on z₁ and z₂
    @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
    nlsolver.z = z₃

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if cache isa Kvaerno3Cache
        @.. broadcast=false z₄=a31 * z₁ + a32 * z₂ + γ * z₃ # use yhat as prediction
    elseif cache isa KenCarp3Cache
        @unpack α41, α42 = cache.tab
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
    end
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = 1
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₄

    step_limiter!(u, integrator, p, t + dt)
    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₄ / dt
end

@muladd function perform_step!(integrator, cache::KenCarp3ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31, α32, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4, ebtilde1, ebtilde2, ebtilde3, ebtilde4 = cache.tab
    alg = unwrap_alg(integrator, true)

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
        z₁ = dt * f(uprev, p, t)
    else
        # FSAL Step 1
        z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation for guess
    nlsolver.z = z₂ = z₁

    nlsolver.tmp = uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * integrator.fsalfirst - z₁
        nlsolver.tmp += ea21 * k1
    end

    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γdt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        z₃ = α31 * z₁ + α32 * z₂
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃
    nlsolver.tmp = tmp
    nlsolver.c = c3

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ = z₂
        u = nlsolver.tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3
    else
        @unpack α41, α42 = cache.tab
        z₄ = α41 * z₁ + α42 * z₂
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

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                  ebtilde1 * k1 + ebtilde2 * k2 + ebtilde3 * k3 + ebtilde4 * k4
        else
            tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

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

@muladd function perform_step!(integrator, cache::KenCarp3Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, k1, k2, k3, k4, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31, α32 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4 = cache.tab
    @unpack ebtilde1, ebtilde2, ebtilde3, ebtilde4 = cache.tab
    alg = unwrap_alg(integrator, true)

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
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation for guess
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + 2γdt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₂
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
                                ea42 * k2 + ea43 * k3
    else
        @unpack α41, α42 = cache.tab
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = 1
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₄
    if integrator.f isa SplitFunction
        f2(k4, u, p, t + dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false u=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 +
                              eb2 * k2 + eb3 * k3 + eb4 * k4
    end

    step_limiter!(u, integrator, p, t + dt)

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ +
                                    btilde4 * z₄ + ebtilde1 * k1 + ebtilde2 * k2 +
                                    ebtilde3 * k3 + ebtilde4 * k4
        else
            @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ +
                                    btilde4 * z₄
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₄ / dt
    end
end

@muladd function perform_step!(integrator, cache::CFNLIRK3ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, c2, c3, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4 = cache.tab
    alg = unwrap_alg(integrator, true)

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
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, k1, k2, k3, k4, atmp, nlsolver = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, c2, c3 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4 = cache.tab

    alg = unwrap_alg(integrator, true)

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
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation for guess
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + c2 * dt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        @.. broadcast=false z₃=z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₂
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
                                ea42 * k2 + ea43 * k3
    else
        @.. broadcast=false z₄=z₂
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = 1
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₄
    if integrator.f isa SplitFunction
        f2(k4, u, p, t + dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false u=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄ + eb1 * k1 +
                              eb2 * k2 + eb3 * k3 + eb4 * k4
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₄ / dt
    end
end

@muladd function perform_step!(integrator, cache::Kvaerno4ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, c3, c4 = cache.tab
    @unpack α21, α31, α32, α41, α42 = cache.tab
    @unpack btilde1, btilde2, btilde3, btilde4, btilde5 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt

    # calculate W
    markfirststage!(nlsolver)

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = zero(u)

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = α31 * z₁ + α32 * z₂

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = α41 * z₁ + α42 * z₂

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use yhat2 for prediction
    nlsolver.z = z₅ = a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₄

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₅ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Kvaerno4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, c3, c4 = cache.tab
    @unpack α21, α31, α32, α41, α42 = cache.tab
    @unpack btilde1, btilde2, btilde3, btilde4, btilde5 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Allow other choices here
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
    nlsolver.z = z₃

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use yhat prediction
    @.. broadcast=false z₅=a41 * z₁ + a42 * z₂ + a43 * z₃ + γ * z₄
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₅

    step_limiter!(u, integrator, p, t + dt)

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₅ / dt
end

@muladd function perform_step!(integrator, cache::KenCarp4ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, c3, c4, c5 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α53, α54, α61, α62, α63, α64, α65 = cache.tab
    @unpack btilde1, btilde3, btilde4, btilde5, btilde6 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65 = cache.tab
    @unpack eb1, eb3, eb4, eb5, eb6 = cache.tab
    @unpack ebtilde1, ebtilde3, ebtilde4, ebtilde5, ebtilde6 = cache.tab
    alg = unwrap_alg(integrator, true)

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

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = z₁

    tmp = uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * integrator.fsalfirst - z₁
        tmp += ea21 * k1
    end
    nlsolver.tmp = tmp
    nlsolver.c = 2γ

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γdt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        z₃ = α31 * z₁ + α32 * z₂
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃
    nlsolver.tmp = tmp
    nlsolver.c = c3

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ = z₂
        u = nlsolver.tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3
    else
        z₄ = α41 * z₁ + α42 * z₂
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄
    nlsolver.tmp = tmp
    nlsolver.c = c4

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ = z₄
        u = nlsolver.tmp + γ * z₄
        k4 = dt * f2(u, p, t + c4 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ + ea51 * k1 + ea52 * k2 +
              ea53 * k3 + ea54 * k4
    else
        z₅ = α51 * z₁ + α52 * z₂ + α53 * z₃ + α54 * z₄
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅
    nlsolver.tmp = tmp
    nlsolver.c = c5

    u = nlsolver.tmp + γ * z₅
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ = z₅
        u = nlsolver.tmp + γ * z₅
        k5 = dt * f2(u, p, t + c5 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + ea61 * k1 + ea62 * k2 +
              ea63 * k3 + ea64 * k4 + ea65 * k5
    else
        z₆ = α61 * z₁ + α62 * z₂ + α63 * z₃ + α64 * z₄ + α65 * z₅
        tmp = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆
    nlsolver.tmp = tmp
    nlsolver.c = 1

    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₆
    if integrator.f isa SplitFunction
        k6 = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + γ * z₆ + eb1 * k1 +
            eb3 * k3 + eb4 * k4 + eb5 * k5 + eb6 * k6
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            tmp = btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ +
                  ebtilde1 * k1 + ebtilde3 * k3 + ebtilde4 * k4 + ebtilde5 * k5 +
                  ebtilde6 * k6
        else
            tmp = btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z₆ ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KenCarp4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, z₆, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack k1, k2, k3, k4, k5, k6 = cache
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, c3, c4, c5 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α53, α54, α61, α62, α63, α64, α65 = cache.tab
    @unpack btilde1, btilde3, btilde4, btilde5, btilde6 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65 = cache.tab
    @unpack eb1, eb3, eb4, eb5, eb6 = cache.tab
    @unpack ebtilde1, ebtilde3, ebtilde4, ebtilde5, ebtilde6 = cache.tab
    alg = unwrap_alg(integrator, true)

    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Allow other choices here
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = 2γ
    markfirststage!(nlsolver)
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + 2γdt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₂
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
                                ea42 * k2 + ea43 * k3
    else
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ .= z₄
        @.. broadcast=false u=tmp + γ * z₄
        f2(k4, u, p, t + c4 * dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ +
                                ea51 * k1 + ea52 * k2 + ea53 * k3 + ea54 * k4
    else
        @.. broadcast=false z₅=α51 * z₁ + α52 * z₂ + α53 * z₃ + α54 * z₄
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅

    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ .= z₅
        @.. broadcast=false u=tmp + γ * z₅
        f2(k5, u, p, t + c5 * dt)
        k5 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ +
                                ea61 * k1 + ea62 * k2 + ea63 * k3 + ea64 * k4 + ea65 * k5
    else
        @.. broadcast=false z₆=α61 * z₁ + α62 * z₂ + α63 * z₃ + α64 * z₄ + α65 * z₅
        @.. broadcast=false tmp=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆

    nlsolver.c = 1
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₆
    if integrator.f isa SplitFunction
        f2(k6, u, p, t + dt)
        k6 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false u=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + γ * z₆ +
                              eb1 * k1 + eb3 * k3 + eb4 * k4 + eb5 * k5 + eb6 * k6
    end

    step_limiter!(u, integrator, p, t + dt)
    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            @.. broadcast=false tmp=btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ +
                                    btilde5 * z₅ + btilde6 * z₆ + ebtilde1 * k1 +
                                    ebtilde3 * k3 + ebtilde4 * k4 + ebtilde5 * k5 +
                                    ebtilde6 * k6
        else
            @.. broadcast=false tmp=btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ +
                                    btilde5 * z₅ + btilde6 * z₆
        end

        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₆ / dt
    end
end

@muladd function perform_step!(integrator, cache::Kvaerno5ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, c3, c4, c5, c6 = cache.tab
    @unpack btilde1, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    @unpack α31, α32, α41, α42, α43, α51, α52, α53, α61, α62, α63 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt

    # calculate W
    markfirststage!(nlsolver)

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = z₁

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = α31 * z₁ + α32 * z₂

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = α41 * z₁ + α42 * z₂ + α43 * z₃

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = z₅ = α51 * z₁ + α52 * z₂ + α53 * z₃

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = z₆ = α61 * z₁ + α62 * z₂ + α63 * z₃

    nlsolver.tmp = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Prediction from embedding
    nlsolver.z = z₇ = a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + γ * z₆

    nlsolver.tmp = uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = 1
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₇

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ +
              btilde7 * z₇
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₇ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Kvaerno5Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, z₆, z₇, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, c3, c4, c5, c6 = cache.tab
    @unpack btilde1, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    @unpack α31, α32, α41, α42, α43, α51, α52, α53, α61, α62, α63 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Allow other choices here
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
    nlsolver.z = z₃

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=α41 * z₁ + α42 * z₂ + α43 * z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    @.. broadcast=false z₅=α51 * z₁ + α52 * z₂ + α53 * z₃
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    @.. broadcast=false z₆=α61 * z₁ + α62 * z₂ + α63 * z₃
    nlsolver.z = z₆

    @.. broadcast=false tmp=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Prediction is embedded method
    @.. broadcast=false z₇=a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + γ * z₆
    nlsolver.z = z₇

    @.. broadcast=false tmp=uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = 1
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₇

    step_limiter!(u, integrator, p, t + dt)
    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
                                btilde6 * z₆ + btilde7 * z₇
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₇ / dt
end

@muladd function perform_step!(integrator, cache::KenCarp5ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, c3, c4, c5, c6, c7 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α61, α62, α71, α72, α73, α74, α75, α81, α82, α83, α84, α85 = cache.tab
    @unpack btilde1, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea43, ea51, ea53, ea54, ea61, ea63, ea64, ea65 = cache.tab
    @unpack ea71, ea73, ea74, ea75, ea76, ea81, ea83, ea84, ea85, ea86, ea87 = cache.tab
    @unpack eb1, eb4, eb5, eb6, eb7, eb8 = cache.tab
    @unpack ebtilde1, ebtilde4, ebtilde5, ebtilde6, ebtilde7, ebtilde8 = cache.tab
    alg = unwrap_alg(integrator, true)

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

    ##### Step 1

    if integrator.f isa SplitFunction
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        z₁ = dt .* f(uprev, p, t)
    else
        # FSAL Step 1
        z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = z₁

    tmp = uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * integrator.fsalfirst - z₁
        tmp += ea21 * k1
    end
    nlsolver.tmp = tmp
    nlsolver.c = 2γ

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γdt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        z₃ = α31 * z₁ + α32 * z₂
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃
    nlsolver.c = c3
    nlsolver.tmp = tmp

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ = z₂
        u = nlsolver.tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a41 * z₁ + a43 * z₃ + ea41 * k1 + ea43 * k3
    else
        z₄ = α41 * z₁ + α42 * z₂
        tmp = uprev + a41 * z₁ + a43 * z₃
    end
    nlsolver.z = z₄
    nlsolver.c = c4
    nlsolver.tmp = tmp

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ = z₂
        u = nlsolver.tmp + γ * z₄
        k4 = dt * f2(u, p, t + c4 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a51 * z₁ + a53 * z₃ + a54 * z₄ + ea51 * k1 + ea53 * k3 + ea54 * k4
    else
        z₅ = α51 * z₁ + α52 * z₂
        tmp = uprev + a51 * z₁ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅
    nlsolver.c = c5
    nlsolver.tmp = tmp

    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ = z₃
        u = nlsolver.tmp + γ * z₅
        k5 = dt * f2(u, p, t + c5 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ + ea61 * k1 + ea63 * k3 +
              ea64 * k4 + ea65 * k5
    else
        z₆ = α61 * z₁ + α62 * z₂
        tmp = uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆
    nlsolver.c = c6
    nlsolver.tmp = tmp

    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ = z₂
        u = nlsolver.tmp + γ * z₆
        k6 = dt * f2(u, p, t + c6 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ + ea71 * k1 +
              ea73 * k3 + ea74 * k4 + ea75 * k5 + ea76 * k6
    else
        z₇ = α71 * z₁ + α72 * z₂ + α73 * z₃ + α74 * z₄ + α75 * z₅
        tmp = uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇
    nlsolver.c = c7
    nlsolver.tmp = tmp

    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    if integrator.f isa SplitFunction
        z₈ = z₅
        u = nlsolver.tmp + γ * z₇
        k7 = dt * f2(u, p, t + c7 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ + ea81 * k1 +
              ea83 * k3 + ea84 * k4 + ea85 * k5 + ea86 * k6 + ea87 * k7
    else
        z₈ = α81 * z₁ + α82 * z₂ + α83 * z₃ + α84 * z₄ + α85 * z₅
        tmp = uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇
    end
    nlsolver.z = z₈
    nlsolver.c = 1
    nlsolver.tmp = tmp

    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₈
    if integrator.f isa SplitFunction
        k8 = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ + γ * z₈ +
            eb1 * k1 + eb4 * k4 + eb5 * k5 + eb6 * k6 + eb7 * k7 + eb8 * k8
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            tmp = btilde1 * z₁ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ +
                  btilde8 * z₈ + ebtilde1 * k1 + ebtilde4 * k4 + ebtilde5 * k5 +
                  ebtilde6 * k6 + ebtilde7 * k7 + ebtilde8 * k8
        else
            tmp = btilde1 * z₁ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ +
                  btilde8 * z₈
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z₈ ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KenCarp5Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈, atmp, nlsolver, step_limiter! = cache
    @unpack k1, k2, k3, k4, k5, k6, k7, k8 = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, c3, c4, c5, c6, c7 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α61, α62, α71, α72, α73, α74, α75, α81, α82, α83, α84, α85 = cache.tab
    @unpack btilde1, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea43, ea51, ea53, ea54, ea61, ea63, ea64, ea65 = cache.tab
    @unpack ea71, ea73, ea74, ea75, ea76, ea81, ea83, ea84, ea85, ea86, ea87 = cache.tab
    @unpack eb1, eb4, eb5, eb6, eb7, eb8 = cache.tab
    @unpack ebtilde1, ebtilde4, ebtilde5, ebtilde6, ebtilde7, ebtilde8 = cache.tab
    alg = unwrap_alg(integrator, true)

    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Allow other choices here
    copyto!(z₂, z₁)
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + 2γdt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        @.. broadcast=false z₃=a31 * z₁ + α32 * z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₃
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a43 * z₃ + ea41 * k1 + ea43 * k3
    else
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
        @.. broadcast=false tmp=uprev + a41 * z₁ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ .= z₂
        @.. broadcast=false u=tmp + γ * z₄
        f2(k4, u, p, t + c4 * dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a51 * z₁ + a53 * z₃ + a54 * z₄ + ea51 * k1 +
                                ea53 * k3 + ea54 * k4
    else
        @.. broadcast=false z₅=α51 * z₁ + α52 * z₂
        @.. broadcast=false tmp=uprev + a51 * z₁ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅

    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ .= z₃
        @.. broadcast=false u=tmp + γ * z₅
        f2(k5, u, p, t + c5 * dt)
        k5 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅ +
                                ea61 * k1 + ea63 * k3 + ea64 * k4 + ea65 * k5
    else
        @.. broadcast=false z₆=α61 * z₁ + α62 * z₂
        @.. broadcast=false tmp=uprev + a61 * z₁ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆

    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ .= z₂
        @.. broadcast=false u=tmp + γ * z₆
        f2(k6, u, p, t + c6 * dt)
        k6 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                                a76 * z₆ + ea71 * k1 + ea73 * k3 + ea74 * k4 + ea75 * k5 +
                                ea76 * k6
    else
        @.. broadcast=false z₇=α71 * z₁ + α72 * z₂ + α73 * z₃ + α74 * z₄ + α75 * z₅
        @.. broadcast=false tmp=uprev + a71 * z₁ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇

    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    if integrator.f isa SplitFunction
        z₈ .= z₅
        @.. broadcast=false u=tmp + γ * z₇
        f2(k7, u, p, t + c7 * dt)
        k7 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ +
                                a87 * z₇ + ea81 * k1 + ea83 * k3 + ea84 * k4 + ea85 * k5 +
                                ea86 * k6 + ea87 * k7
    else
        @.. broadcast=false z₈=α81 * z₁ + α82 * z₂ + α83 * z₃ + α84 * z₄ + α85 * z₅
        @.. broadcast=false tmp=uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇
    end
    nlsolver.z = z₈

    nlsolver.c = 1
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₈
    if integrator.f isa SplitFunction
        f2(k8, u, p, t + dt)
        k8 .*= dt
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false u=uprev + a81 * z₁ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ +
                              γ * z₈ + eb1 * k1 + eb4 * k4 + eb5 * k5 + eb6 * k6 +
                              eb7 * k7 + eb8 * k8
    end

    step_limiter!(u, integrator, p, t + dt)
    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            @.. broadcast=false tmp=btilde1 * z₁ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈ +
                                    ebtilde1 * k1 + ebtilde4 * k4 + ebtilde5 * k5 +
                                    ebtilde6 * k6 + ebtilde7 * k7 + ebtilde8 * k8
        else
            @.. broadcast=false tmp=btilde1 * z₁ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈
        end

        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₈ / dt
    end
end

@muladd function perform_step!(integrator, cache::KenCarp47ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a73, a74, a75, a76, c3, c4, c5, c6 = cache.tab
    @unpack α31, α32, α41, α42, α43, α51, α52, α61, α62, α63, α71, α72, α73, α74, α75, α76 = cache.tab
    @unpack btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65, ea71, ea72, ea73, ea74, ea75, ea76 = cache.tab
    @unpack eb3, eb4, eb5, eb6, eb7 = cache.tab
    @unpack ebtilde3, ebtilde4, ebtilde5, ebtilde6, ebtilde7 = cache.tab
    alg = unwrap_alg(integrator, true)

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

    ##### Step 1

    if integrator.f isa SplitFunction
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        z₁ = dt .* f(uprev, p, t)
    else
        # FSAL Step 1
        z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = z₁

    tmp = uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * integrator.fsalfirst - z₁
        tmp += ea21 * k1
    end
    nlsolver.tmp = tmp
    nlsolver.c = 2γ

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γdt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        z₃ = α31 * z₁ + α32 * z₂
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
        z₄ = α41 * z₁ + α42 * z₂ + α43 * z₃
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄
    nlsolver.tmp = tmp
    nlsolver.c = c4

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ = z₁
        u = nlsolver.tmp + γ * z₄
        k4 = dt * f2(u, p, t + c4 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ + ea51 * k1 + ea52 * k2 +
              ea53 * k3 + ea54 * k4
    else
        z₅ = α51 * z₁ + α52 * z₂
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅
    nlsolver.tmp = tmp
    nlsolver.c = c5

    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ = z₃
        u = nlsolver.tmp + γ * z₅
        k5 = dt * f2(u, p, t + c5 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅ + ea61 * k1 +
              ea62 * k2 + ea63 * k3 + ea64 * k4 + ea65 * k5
    else
        z₆ = α61 * z₁ + α62 * z₂ + α63 * z₃
        tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆
    nlsolver.tmp = tmp
    nlsolver.c = c6

    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ = z₆
        u = nlsolver.tmp + γ * z₆
        k6 = dt * f2(u, p, t + c6 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ + ea71 * k1 + ea72 * k2 +
              ea73 * k3 + ea74 * k4 + ea75 * k5 + ea76 * k6
    else
        z₇ = α71 * z₁ + α72 * z₂ + α73 * z₃ + α74 * z₄ + α75 * z₅ + +α76 * z₆
        tmp = uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇
    nlsolver.c = 1
    nlsolver.tmp = tmp

    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₇
    if integrator.f isa SplitFunction
        k7 = dt * f2(u, p, t + dt)
        u = uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ + γ * z₇ + eb3 * k3 +
            eb4 * k4 + eb5 * k5 + eb6 * k6 + eb7 * k7
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            tmp = btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ +
                  ebtilde3 * k3 + ebtilde4 * k4 + ebtilde5 * k5 + ebtilde6 * k6 +
                  ebtilde7 * k7
        else
            tmp = btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z₇ ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KenCarp47Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, z₆, z₇, atmp, nlsolver = cache
    @unpack k1, k2, k3, k4, k5, k6, k7 = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a73, a74, a75, a76, c3, c4, c5, c6 = cache.tab
    @unpack α31, α32, α41, α42, α43, α51, α52, α61, α62, α63, α71, α72, α73, α74, α75, α76 = cache.tab
    @unpack btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65, ea71, ea72, ea73, ea74, ea75, ea76 = cache.tab
    @unpack eb3, eb4, eb5, eb6, eb7 = cache.tab
    @unpack ebtilde3, ebtilde4, ebtilde5, ebtilde6, ebtilde7 = cache.tab
    alg = unwrap_alg(integrator, true)

    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Allow other choices here
    z₂ .= z₁
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + 2γdt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        #Guess is from Hermite derivative on z₁ and z₂
        @.. broadcast=false z₃=a31 * z₁ + α32 * z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₃
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
                                ea42 * k2 + ea43 * k3
    else
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂ + α43 * z₃
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ .= z₁
        @.. broadcast=false u=tmp + γ * z₄
        f2(k4, u, p, t + c4 * dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ +
                                ea51 * k1 + ea52 * k2 + ea53 * k3 + ea54 * k4
    else
        @.. broadcast=false z₅=α51 * z₁ + α52 * z₂
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅

    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ .= z₃
        @.. broadcast=false u=tmp + γ * z₅
        f2(k5, u, p, t + c5 * dt)
        k5 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ +
                                a65 * z₅ + ea61 * k1 + ea62 * k2 + ea63 * k3 + ea64 * k4 +
                                ea65 * k5
    else
        @.. broadcast=false z₆=α61 * z₁ + α62 * z₂ + α63 * z₃
        @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆

    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ .= z₆
        @.. broadcast=false u=tmp + γ * z₆
        f2(k6, u, p, t + c6 * dt)
        k6 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ +
                                ea71 * k1 + ea72 * k2 + ea73 * k3 + ea74 * k4 + ea75 * k5 +
                                ea76 * k6
    else
        @.. broadcast=false z₇=α71 * z₁ + α72 * z₂ + α73 * z₃ + α74 * z₄ + α75 * z₅ +
                               α76 * z₆
        @.. broadcast=false tmp=uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇

    nlsolver.c = 1
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₇
    if integrator.f isa SplitFunction
        f2(k7, u, p, t + dt)
        k7 .*= dt
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false u=uprev + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ + γ * z₇ +
                              eb3 * k3 + eb4 * k4 + eb5 * k5 + eb6 * k6 + eb7 * k7
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            @.. broadcast=false tmp=btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇ + ebtilde3 * k3 +
                                    ebtilde4 * k4 + ebtilde5 * k5 + ebtilde6 * k6 +
                                    ebtilde7 * k7
        else
            @.. broadcast=false tmp=btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇
        end

        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₇ / dt
    end
end

@muladd function perform_step!(integrator, cache::KenCarp58ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    nlsolver = cache.nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a83, a84, a85, a86, a87, c3, c4, c5, c6, c7 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α61, α62, α63, α71, α72, α73, α81, α82, α83, α84, α85, α86, α87 = cache.tab
    @unpack btilde3, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65 = cache.tab
    @unpack ea71, ea72, ea73, ea74, ea75, ea76, ea81, ea82, ea83, ea84, ea85, ea86, ea87 = cache.tab
    @unpack eb3, eb4, eb5, eb6, eb7, eb8 = cache.tab
    @unpack ebtilde3, ebtilde4, ebtilde5, ebtilde6, ebtilde7, ebtilde8 = cache.tab
    alg = unwrap_alg(integrator, true)

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

    ##### Step 1

    if integrator.f isa SplitFunction
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        z₁ = dt .* f(uprev, p, t)
    else
        # FSAL Step 1
        z₁ = dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Add extrapolation choice

    nlsolver.z = z₂ = z₁

    tmp = uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        k1 = dt * integrator.fsalfirst - z₁
        tmp += ea21 * k1
    end
    nlsolver.tmp = tmp
    nlsolver.c = 2γ

    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ = z₂
        u = nlsolver.tmp + γ * z₂
        k2 = dt * f2(u, p, t + 2γdt)
        integrator.stats.nf2 += 1
        tmp = uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        z₃ = α31 * z₁ + α32 * z₂
        tmp = uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃
    nlsolver.c = c3
    nlsolver.tmp = tmp

    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ = z₁
        u = nlsolver.tmp + γ * z₃
        k3 = dt * f2(u, p, t + c3 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 + ea42 * k2 + ea43 * k3
    else
        z₄ = α41 * z₁ + α42 * z₂
        tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄
    nlsolver.c = c4
    nlsolver.tmp = tmp

    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ = z₂
        u = nlsolver.tmp + γ * z₄
        k4 = dt * f2(u, p, t + c4 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ + ea51 * k1 + ea52 * k2 +
              ea53 * k3 + ea54 * k4
    else
        z₅ = α51 * z₁ + α52 * z₂
        tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅
    nlsolver.c = c5
    nlsolver.tmp = tmp

    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ = z₃
        u = nlsolver.tmp + γ * z₅
        k5 = dt * f2(u, p, t + c5 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅ + ea61 * k1 +
              ea62 * k2 + ea63 * k3 + ea64 * k4 + ea65 * k5
    else
        z₆ = α61 * z₁ + α62 * z₂ + α63 * z₃
        tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆
    nlsolver.c = c6
    nlsolver.tmp = tmp

    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ = z₃
        u = nlsolver.tmp + γ * z₆
        k6 = dt * f2(u, p, t + c6 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆ +
              ea71 * k1 + ea72 * k2 + ea73 * k3 + ea74 * k4 + ea75 * k5 + ea76 * k6
    else
        z₇ = α71 * z₁ + α72 * z₂ + α73 * z₃
        tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇
    nlsolver.c = c7
    nlsolver.tmp = tmp

    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    if integrator.f isa SplitFunction
        z₈ = z₇
        u = nlsolver.tmp + γ * z₇
        k7 = dt * f2(u, p, t + c7 * dt)
        integrator.stats.nf2 += 1
        tmp = uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ + ea81 * k1 +
              ea82 * k2 + ea83 * k3 + ea84 * k4 + ea85 * k5 + ea86 * k6 + ea87 * k7
    else
        z₈ = α81 * z₁ + α82 * z₂ + α83 * z₃ + α84 * z₄ + α85 * z₅ + α86 * z₆ + α87 * z₇
        tmp = uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇
    end
    nlsolver.z = z₈
    nlsolver.c = 1
    nlsolver.tmp = tmp

    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₈
    if integrator.f isa SplitFunction
        k8 = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ + γ * z₈ +
            eb3 * k3 + eb4 * k4 + eb5 * k5 + eb6 * k6 + eb7 * k7 + eb8 * k8
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            tmp = btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ +
                  btilde8 * z₈ + ebtilde3 * k3 + ebtilde4 * k4 + ebtilde5 * k5 +
                  ebtilde6 * k6 + ebtilde7 * k7 + ebtilde8 * k8
        else
            tmp = btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ +
                  btilde8 * z₈
        end
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z₈ ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KenCarp58Cache, repeat_step = false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈, atmp, nlsolver = cache
    @unpack k1, k2, k3, k4, k5, k6, k7, k8 = cache
    @unpack tmp = nlsolver
    @unpack γ, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a83, a84, a85, a86, a87, c3, c4, c5, c6, c7 = cache.tab
    @unpack α31, α32, α41, α42, α51, α52, α61, α62, α63, α71, α72, α73, α81, α82, α83, α84, α85, α86, α87 = cache.tab
    @unpack btilde3, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    @unpack ea21, ea31, ea32, ea41, ea42, ea43, ea51, ea52, ea53, ea54, ea61, ea62, ea63, ea64, ea65 = cache.tab
    @unpack ea71, ea72, ea73, ea74, ea75, ea76, ea81, ea82, ea83, ea84, ea85, ea86, ea87 = cache.tab
    @unpack eb3, eb4, eb5, eb6, eb7, eb8 = cache.tab
    @unpack ebtilde3, ebtilde4, ebtilde5, ebtilde6, ebtilde7, ebtilde8 = cache.tab
    alg = unwrap_alg(integrator, true)

    if integrator.f isa SplitFunction
        f = integrator.f.f1
        f2 = integrator.f.f2
    else
        f = integrator.f
    end

    # precalculations
    γdt = γ * dt

    markfirststage!(nlsolver)

    ##### Step 1

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        # Explicit tableau is not FSAL
        # Make this not compute on repeat
        f(z₁, integrator.uprev, p, integrator.t)
        z₁ .*= dt
    else
        # FSAL Step 1
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    end

    ##### Step 2

    # TODO: Allow other choices here
    z₂ .= z₁
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁

    if integrator.f isa SplitFunction
        # This assumes the implicit part is cheaper than the explicit part
        @.. broadcast=false k1=dt * integrator.fsalfirst - z₁
        @.. broadcast=false tmp+=ea21 * k1
    end

    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    if integrator.f isa SplitFunction
        z₃ .= z₂
        @.. broadcast=false u=tmp + γ * z₂
        f2(k2, u, p, t + 2γdt)
        k2 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂ + ea31 * k1 + ea32 * k2
    else
        # Guess is from Hermite derivative on z₁ and z₂
        @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
        @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    end
    nlsolver.z = z₃

    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    if integrator.f isa SplitFunction
        z₄ .= z₁
        @.. broadcast=false u=tmp + γ * z₃
        f2(k3, u, p, t + c3 * dt)
        k3 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃ + ea41 * k1 +
                                ea42 * k2 + ea43 * k3
    else
        @.. broadcast=false z₄=α41 * z₁ + α42 * z₂
        @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    end
    nlsolver.z = z₄

    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    if integrator.f isa SplitFunction
        z₅ .= z₂
        @.. broadcast=false u=tmp + γ * z₄
        f2(k4, u, p, t + c4 * dt)
        k4 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄ +
                                ea51 * k1 + ea52 * k2 + ea53 * k3 + ea54 * k4
    else
        @.. broadcast=false z₅=α51 * z₁ + α52 * z₂
        @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    end
    nlsolver.z = z₅

    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    if integrator.f isa SplitFunction
        z₆ .= z₃
        @.. broadcast=false u=tmp + γ * z₅
        f2(k5, u, p, t + c5 * dt)
        k5 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ +
                                a65 * z₅ + ea61 * k1 + ea62 * k2 + ea63 * k3 + ea64 * k4 +
                                ea65 * k5
    else
        @.. broadcast=false z₆=α61 * z₁ + α62 * z₂ + α63 * z₃
        @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    end
    nlsolver.z = z₆

    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    if integrator.f isa SplitFunction
        z₇ .= z₃
        @.. broadcast=false u=tmp + γ * z₆
        f2(k6, u, p, t + c6 * dt)
        k6 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ +
                                a75 * z₅ + a76 * z₆ + ea71 * k1 + ea72 * k2 + ea73 * k3 +
                                ea74 * k4 + ea75 * k5 + ea76 * k6
    else
        @.. broadcast=false z₇=α71 * z₁ + α72 * z₂ + α73 * z₃
        @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ +
                                a75 * z₅ + a76 * z₆
    end
    nlsolver.z = z₇

    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    if integrator.f isa SplitFunction
        z₈ .= z₇
        @.. broadcast=false u=tmp + γ * z₇
        f2(k7, u, p, t + c7 * dt)
        k7 .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false tmp=uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ +
                                a87 * z₇ + ea81 * k1 + ea82 * k2 + ea83 * k3 + ea84 * k4 +
                                ea85 * k5 + ea86 * k6 + ea87 * k7
    else
        @.. broadcast=false z₈=α81 * z₁ + α82 * z₂ + α83 * z₃ + α84 * z₄ + α85 * z₅ +
                               α86 * z₆ + α87 * z₇
        @.. broadcast=false tmp=uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇
    end
    nlsolver.z = z₈

    nlsolver.c = 1
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₈
    if integrator.f isa SplitFunction
        f2(k8, u, p, t + dt)
        k8 .*= dt
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false u=uprev + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇ +
                              γ * z₈ + eb3 * k3 + eb4 * k4 + eb5 * k5 + eb6 * k6 +
                              eb7 * k7 + eb8 * k8
    end

    ################################### Finalize

    if integrator.opts.adaptive
        if integrator.f isa SplitFunction
            @.. broadcast=false tmp=btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈ +
                                    ebtilde3 * k3 + ebtilde4 * k4 + ebtilde5 * k5 +
                                    ebtilde6 * k6 + ebtilde7 * k7 + ebtilde8 * k8
        else
            @.. broadcast=false tmp=btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
                                    btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈
        end

        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))

            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast=z₈ / dt
    end
end

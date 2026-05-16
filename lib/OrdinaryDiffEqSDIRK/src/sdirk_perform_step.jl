function initialize!(integrator, cache::SDIRKConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::SDIRKMutableCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(
        integrator, cache::ImplicitEulerConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else # :constant
        nlsolver.z = zero(u)
    end

    nlsolver.tmp = uprev
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + z

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        tmp = r *
            integrator.opts.internalnorm.((u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    else
        OrdinaryDiffEqCore.set_EEst!(integrator, 1)
    end

    integrator.fsallast = f(u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
            integrator.opts.abstol
        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + (integrator.opts.internalnorm(atmp, t)))
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; atmp, nlsolver, step_limiter!) = cache
    (; z, tmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast = false z = dt * integrator.fsalfirst
    else # :constant
        z .= zero(eltype(u))
    end

    nlsolver.tmp .= uprev
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = uprev + z

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        @.. broadcast = false tmp = r * integrator.opts.internalnorm(
            (u - uprev) / dt1 -
                (uprev - uprev2) / dt2, t
        )
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    else
        OrdinaryDiffEqCore.set_EEst!(integrator, 1)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(integrator.fsallast, u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        @.. broadcast = false atmp = ifelse(cache.algebraic_vars, integrator.fsallast, false) /
            integrator.opts.abstol
        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + (integrator.opts.internalnorm(atmp, t)))
    end
end

@muladd function perform_step!(
        integrator, cache::TrapezoidConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    # precalculations
    γ = 1 // 2
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess: constant extrapolation
    nlsolver.z = uprev

    if f.mass_matrix === I
        nlsolver.tmp = @.. broadcast = false uprev * inv(γdt) + integrator.fsalfirst
    else
        nlsolver.tmp = (f.mass_matrix * uprev) .* inv(γdt) .+ integrator.fsalfirst
    end
    nlsolver.α = 1
    nlsolver.γ = γ
    nlsolver.method = COEFFICIENT_MULTISTEP
    u = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            # tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
            DD31 = (u - uprev) / dt1 - (uprev - uprev2) / dt2
            DD30 = (uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4
            tmp = r * integrator.opts.internalnorm((DD31 - DD30) / dt5, t)
            atmp = calculate_residuals(
                tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
            if OrdinaryDiffEqCore.get_EEst(integrator) <= 1
                cache.uprev3 = uprev2
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
            cache.uprev3 = integrator.uprev2
            cache.tprev2 = integrator.tprev
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; atmp, nlsolver, step_limiter!) = cache
    (; z, tmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    γ = 1 // 2
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess: constant extrapolation
    @.. broadcast = false z = uprev
    invγdt = inv(γdt)
    if mass_matrix === I
        @.. broadcast = false tmp = uprev * invγdt + integrator.fsalfirst
    else
        mul!(u, mass_matrix, uprev)
        @.. broadcast = false tmp = u * invγdt + integrator.fsalfirst
    end
    nlsolver.α = 1
    nlsolver.γ = γ
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            # @.. broadcast=false tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
            @.. broadcast = false tmp = r * integrator.opts.internalnorm(
                (
                    (
                        (u - uprev) / dt1 -
                            (uprev - uprev2) / dt2
                    ) #DD31
                        -
                        (
                        (uprev - uprev2) / dt3 -
                            (uprev2 - uprev3) /
                            dt4
                    )
                ) /
                    dt5,
                t
            )
            calculate_residuals!(
                atmp, tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
            if OrdinaryDiffEqCore.get_EEst(integrator) <= 1
                copyto!(cache.uprev3, uprev2)
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
            copyto!(cache.uprev3, integrator.uprev2)
            cache.tprev2 = integrator.tprev
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(integrator.fsallast, u, p, t + dt)
end

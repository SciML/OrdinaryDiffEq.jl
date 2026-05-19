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

function calculate_error_estimate!(
        integrator, cache::ImplicitEulerConstantCache, ::Nothing, t
    )
    (; uprev, u, dt) = integrator
    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
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
    integrator.u = u

    calculate_error_estimate!(integrator, cache, nothing, t)

    integrator.fsallast = f(u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
            integrator.opts.abstol
        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + (integrator.opts.internalnorm(atmp, t)))
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function calculate_error_estimate!(
        integrator, cache::ImplicitEulerCache, ::Nothing, t
    )
    (; uprev, u, dt) = integrator
    (; atmp, nlsolver) = cache
    (; tmp) = nlsolver
    if integrator.opts.adaptive && integrator.success_iter > 0
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
        @.. broadcast = false tmp = r * integrator.opts.internalnorm(
            (u - uprev) / dt1 - (uprev - uprev2) / dt2, t
        )
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    else
        OrdinaryDiffEqCore.set_EEst!(integrator, 1)
    end
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
    integrator.u = u

    step_limiter!(u, integrator, p, t + dt)

    calculate_error_estimate!(integrator, cache, nothing, t)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(integrator.fsallast, u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        @.. broadcast = false atmp = ifelse(cache.algebraic_vars, integrator.fsallast, false) /
            integrator.opts.abstol
        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + (integrator.opts.internalnorm(atmp, t)))
    end
end

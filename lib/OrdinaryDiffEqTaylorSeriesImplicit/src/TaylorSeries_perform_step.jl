function initialize!(integrator, cache::ImplicitTaylorConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitTaylorMutableCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ImplicitTaylor1ConstantCache,
        repeat_step = false)
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
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end

    integrator.fsallast = f(u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                  integrator.opts.abstol
        integrator.EEst += integrator.opts.internalnorm(atmp, t)
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ImplicitTaylor1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; μ, κ, tmp, atmp, step_limiter!) = cache
    alg = unwrap_alg(integrator, true)
    (; internalnorm, abstol, reltol, adaptive) = integrator.opts
    (; maxiters) = alg

    # the reshape should really be moved to build_jet
    nlf = (u_next_vec) -> begin
        u_next = reshape(u_next_vec, size(uprev))
        utaylor = reshape(cache.jet(u_next, t + μ * dt), size(uprev))
        du = get_coefficient(utaylor, 1)
        _vec(u_next .- uprev .- μ * dt .* du)
    end
    #=cache.tmp .= u
    nlprob = NonlinearProblem(nlf, cache.tmp, u)
    sol = solve(nlprob)=#


    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))
    fail_convergence = true
    iter = 0
    dw = cache.linsolve.u
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1
        # evaluate function
        tmp .= reshape(nlf(_vec(u)), size(tmp))
        #@show tmp
        linsolve = cache.linsolve
        needfactor = iter == 1 # && newW once we aren't using ForwardDiff
        A =  needfactor ? ForwardDiff.jacobian(nlf, _vec(u)) : nothing
        needfactor && @assert A isa Matrix
        linres = dolinsolve(integrator, linsolve; A, b = _vec(tmp),
            linu = _vec(dw))
        #@show dw
        cache.linsolve = linres.cache
        integrator.stats.nsolve += 1

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        calculate_residuals!(atmp, dw, uprev, u, abstol, reltol, internalnorm, t)
        ndw = internalnorm(atmp, t)

        u .-= reshape(dw, size(u))
        # check divergence (not in initial step)
        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 2) && (cache.status = Divergence)
            if diverge
                break
            end
        end
        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = Convergence
            fail_convergence = false
            break
        end
    end
    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    if μ != one(μ)
        utaylor = reshape(cache.jet(u, t + μ * dt), size(u))
        u1 = get_coefficient(utaylor, 1)
        u .+= (1 - μ) * dt * u1
    end
    # step_limiter!(u, integrator, p, t + dt)

    if adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        @.. broadcast=false tmp=r * integrator.opts.internalnorm(
            (u - uprev) / dt1 -
            (uprev - uprev2) / dt2, t)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(integrator.fsallast, u, p, t + dt)
end

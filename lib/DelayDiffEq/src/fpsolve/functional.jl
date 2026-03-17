function OrdinaryDiffEqNonlinearSolve.compute_step!(
        fpsolver::FPSolver{<:NLFunctional},
        integrator::DDEIntegrator
    )
    # update ODE integrator to next time interval together with correct interpolation
    if fpsolver.iter == 1
        advance_ode_integrator!(integrator)
    else
        update_ode_integrator!(integrator)
    end

    return compute_step_fixedpoint!(fpsolver, integrator)
end

function OrdinaryDiffEqNonlinearSolve.compute_step!(
        fpsolver::FPSolver{<:NLAnderson, false},
        integrator::DDEIntegrator
    )
    (; cache, iter) = fpsolver
    (; aa_start) = cache

    # perform Anderson acceleration
    previter = iter - 1
    if previter == aa_start
        # update cached values for next step of Anderson acceleration
        cache.dzold = cache.dz
        cache.z₊old = integrator.u
    elseif previter > aa_start
        # actually perform Anderson acceleration
        integrator.u = OrdinaryDiffEqNonlinearSolve.anderson(integrator.u, cache)
        integrator.stats.nsolve += 1
    end

    # update ODE integrator to next time interval together with correct interpolation
    if iter == 1
        advance_ode_integrator!(integrator)
    else
        # force recomputation of all interpolation data if Anderson acceleration was performed
        update_ode_integrator!(integrator, previter > aa_start)
    end

    # compute next step
    return compute_step_fixedpoint!(fpsolver, integrator)
end

function OrdinaryDiffEqNonlinearSolve.compute_step!(
        fpsolver::FPSolver{<:NLAnderson, true},
        integrator::DDEIntegrator
    )
    (; cache, iter) = fpsolver
    (; aa_start) = cache

    # perform Anderson acceleration
    previter = iter - 1
    if previter == aa_start
        # update cached values for next step of Anderson acceleration
        @.. cache.dzold = cache.dz
        @.. cache.z₊old = integrator.u
    elseif previter > aa_start
        # actually perform Anderson acceleration
        OrdinaryDiffEqNonlinearSolve.anderson!(integrator.u, cache)
        integrator.stats.nsolve += 1
    end

    # update ODE integrator to next time interval together with correct interpolation
    if iter == 1
        advance_ode_integrator!(integrator)
    else
        # force recomputation of all interpolation data if Anderson acceleration was performed
        update_ode_integrator!(integrator, previter > aa_start)
    end

    # compute next step
    return compute_step_fixedpoint!(fpsolver, integrator)
end

function compute_step_fixedpoint!(
        fpsolver::FPSolver{
            <:Union{NLFunctional, NLAnderson},
            false,
        },
        integrator::DDEIntegrator
    )
    (; t, opts) = integrator
    (; cache) = fpsolver
    ode_integrator = integrator.integrator

    # recompute next integration step
    OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache, true)

    # compute residuals
    dz = integrator.u .- ode_integrator.u
    atmp = DiffEqBase.calculate_residuals(
        dz, ode_integrator.u, integrator.u,
        opts.abstol, opts.reltol, opts.internalnorm,
        t
    )

    # cache results
    if isdefined(cache, :dz)
        cache.dz = dz
    end

    return opts.internalnorm(atmp, t)
end

function compute_step_fixedpoint!(
        fpsolver::FPSolver{
            <:Union{NLFunctional, NLAnderson},
            true,
        },
        integrator::DDEIntegrator
    )
    (; t, opts) = integrator
    (; cache) = fpsolver
    (; dz, atmp) = cache
    ode_integrator = integrator.integrator

    # recompute next integration step
    OrdinaryDiffEqCore.perform_step!(integrator, integrator.cache, true)

    # compute residuals
    @.. dz = integrator.u - ode_integrator.u
    DiffEqBase.calculate_residuals!(
        atmp, dz, ode_integrator.u, integrator.u,
        opts.abstol, opts.reltol, opts.internalnorm, t
    )

    residual = opts.internalnorm(atmp, t)
    @SciMLMessage(
        lazy"Fixed-point iteration residual = $residual at t = $t",
        opts.verbose, :residual_control
    )
    return residual
end

## resize!

function Base.resize!(fpcache::FPFunctionalCache, i::Int)
    resize!(fpcache.atmp, i)
    resize!(fpcache.dz, i)
    return nothing
end

function Base.resize!(
        fpcache::FPAndersonCache, fpsolver::FPSolver{<:NLAnderson},
        integrator::DDEIntegrator, i::Int
    )
    return resize!(fpcache, fpsolver.alg, i)
end

function Base.resize!(fpcache::FPAndersonCache, fpalg::NLAnderson, i::Int)
    (; z₊old, Δz₊s) = fpcache

    resize!(fpcache.atmp, i)
    resize!(fpcache.dz, i)
    resize!(fpcache.dzold, i)
    resize!(z₊old, i)

    # update history of Anderson cache
    max_history_old = length(Δz₊s)
    max_history = min(fpalg.max_history, fpalg.max_iter, i)

    resize!(fpcache.γs, max_history)
    resize!(fpcache.Δz₊s, max_history)

    if max_history != max_history_old
        fpcache.Q = typeof(fpcache.Q)(undef, i, max_history)
        fpcache.R = typeof(fpcache.R)(undef, max_history, max_history)
    end

    max_history = length(Δz₊s)
    if max_history > max_history_old
        for i in (max_history_old + 1):max_history
            Δz₊s[i] = zero(z₊old)
        end
    end

    return nothing
end

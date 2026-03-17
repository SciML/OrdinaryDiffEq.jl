@muladd function perform_step!(integrator, cache::TauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    tmp = c(uprev, p, t, P.dW, nothing)
    integrator.u = uprev .+ tmp

    if integrator.opts.adaptive
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            newrate = P.cache.rate(integrator.u, p, t + dt)
            EEstcache = @. abs(newrate - oldrate) /
                max(50integrator.opts.reltol * oldrate, integrator.rate_constants / integrator.dt)
            integrator.EEst = maximum(EEstcache)
            if integrator.EEst <= 1
                P.cache.currate = newrate
            end
        elseif integrator.alg isa CaoTauLeaping
            # Calculate τ as EEst
        end
    end
end

@muladd function perform_step!(integrator, cache::TauLeapingCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; tmp, newrate, EEstcache) = cache
    c(tmp, uprev, p, t, P.dW, nothing)
    @.. u = uprev + tmp

    if integrator.opts.adaptive
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            P.cache.rate(newrate, u, p, t + dt)
            @.. EEstcache = abs(newrate - oldrate) /
                max(50integrator.opts.reltol * oldrate, integrator.rate_constants / integrator.dt)
            integrator.EEst = maximum(EEstcache)
            if integrator.EEst <= 1
                P.cache.currate .= newrate
            end
        elseif integrator.alg isa CaoTauLeaping
            # Calculate τ as EEst
        end
    end
end

# ImplicitTauLeaping: First-order implicit (backward Euler) tau-leaping
# Based on Rathinam, Petzold, Cao, Gillespie (2003)
#
# Implicit equation:
#   X_{n+1} = X_n + ν*k + dt*(drift(X_{n+1}) - drift(X_n))
# where k ~ Poisson(dt * a(X_n)) and drift(u) = ν*a(u) = c(u, p, t, rate(u, p, t), nothing)
#
# Rearranged for nlsolver:
#   X_{n+1} = tmp + z
# where:
#   tmp = X_n + ν*k - dt*drift(X_n)  (explicit contribution)
#   z = dt*drift(X_{n+1})            (solved by nlsolver)
#
# Uses standard nlsolver infrastructure from OrdinaryDiffEqNonlinearSolve

@muladd function perform_step!(integrator, cache::ImplicitTauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; poisson_counts, rate_at_uprev, nlsolver) = cache
    rng = P.rng
    repeat_step = false

    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    # Step 1: Get rates at current state and generate Poisson counts
    rate_at_uprev = P.cache.rate(uprev, p, t)
    poisson_counts = JumpProcesses.pois_rand.(Ref(rng), dt .* rate_at_uprev)

    # Step 2: Compute explicit contributions
    # jump_contribution = ν*k
    jump_contribution = c(uprev, p, t, poisson_counts, nothing)
    # drift_at_uprev = ν*a(X_n)
    drift_at_uprev = c(uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up nlsolver
    # tmp = X_n + ν*k - dt*drift(X_n)
    nlsolver.tmp = uprev .+ jump_contribution .- dt .* drift_at_uprev
    # Initial guess: z = dt * drift(uprev) (linear extrapolation)
    nlsolver.z = dt .* drift_at_uprev

    # Step 4: Solve using nlsolver
    # The nlsolver solves: z = dt * drift(tmp + z)
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    # Step 5: Final update: X_{n+1} = tmp + z
    integrator.u = nlsolver.tmp .+ z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.currate = P.cache.rate(integrator.u, p, t + dt)
    end

    return nothing
end

@muladd function perform_step!(integrator, cache::ImplicitTauLeapingCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; poisson_counts, rate_at_uprev, nlsolver) = cache
    rng = P.rng
    repeat_step = false

    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    # Step 1: Get rates at current state (in-place) and generate Poisson counts
    P.cache.rate(rate_at_uprev, uprev, p, t)
    @. poisson_counts = JumpProcesses.pois_rand(rng, dt * rate_at_uprev)

    # Step 2: Compute explicit contributions
    # Use nlsolver.tmp for intermediate storage
    # First compute jump_contribution = ν*k into nlsolver.tmp
    c(nlsolver.tmp, uprev, p, t, poisson_counts, nothing)

    # Compute drift at uprev: drift(X_n) = ν*a(X_n) into nlsolver.z (will be overwritten)
    c(nlsolver.z, uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up nlsolver
    # tmp = X_n + ν*k - dt*drift(X_n)
    @.. nlsolver.tmp = uprev + nlsolver.tmp - dt * nlsolver.z
    # Initial guess: z = dt * drift(uprev) (linear extrapolation)
    @.. nlsolver.z = dt * nlsolver.z

    # Step 4: Solve using nlsolver
    # The nlsolver solves: z = dt * drift(tmp + z)
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    # Step 5: Final update: X_{n+1} = tmp + z
    @.. u = nlsolver.tmp + z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.rate(P.cache.currate, u, p, t + dt)
    end

    return nothing
end

# ThetaTrapezoidalTauLeaping: Implicit weak second order tau-leaping
# Based on Hu, Li, Min (2011) and Anderson, Mattingly (2011)
#
# Implicit equation:
#   X_{n+1} = X_n + ν*k + θ*dt*(drift(X_{n+1}) - drift(X_n))
# where k ~ Poisson(dt * a(X_n)) and drift(u) = ν*a(u) = c(u, p, t, rate(u, p, t), nothing)
#
# Rearranged for nlsolver:
#   X_{n+1} = tmp + θ*z
# where:
#   tmp = X_n + ν*k - θ*dt*drift(X_n)  (explicit contribution)
#   z = dt*drift(X_{n+1})              (solved by nlsolver)
#
# Uses standard nlsolver infrastructure from OrdinaryDiffEqNonlinearSolve

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; poisson_counts, rate_at_uprev, theta, nlsolver) = cache
    rng = P.rng
    repeat_step = false

    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    # Step 1: Get rates at current state and generate Poisson counts
    rate_at_uprev = P.cache.rate(uprev, p, t)
    poisson_counts = JumpProcesses.pois_rand.(Ref(rng), dt .* rate_at_uprev)

    # Step 2: Compute explicit contributions
    # jump_contribution = ν*k
    jump_contribution = c(uprev, p, t, poisson_counts, nothing)
    # drift_at_uprev = ν*a(X_n)
    drift_at_uprev = c(uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up nlsolver
    # tmp = X_n + ν*k - θ*dt*drift(X_n)
    nlsolver.tmp = uprev .+ jump_contribution .- theta .* dt .* drift_at_uprev
    # Initial guess: z = dt * drift(uprev) (linear extrapolation)
    nlsolver.z = dt .* drift_at_uprev

    # Step 4: Solve using nlsolver
    # The nlsolver solves: z = dt * drift(tmp + θ*z)
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    # Step 5: Final update: X_{n+1} = tmp + θ*z
    integrator.u = nlsolver.tmp .+ theta .* z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.currate = P.cache.rate(integrator.u, p, t + dt)
    end

    return nothing
end

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; poisson_counts, rate_at_uprev, theta, nlsolver) = cache
    rng = P.rng
    repeat_step = false

    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    # Step 1: Get rates at current state (in-place) and generate Poisson counts
    P.cache.rate(rate_at_uprev, uprev, p, t)
    @. poisson_counts = JumpProcesses.pois_rand(rng, dt * rate_at_uprev)

    # Step 2: Compute explicit contributions
    # Use nlsolver.tmp for intermediate storage
    # First compute jump_contribution = ν*k into nlsolver.tmp
    c(nlsolver.tmp, uprev, p, t, poisson_counts, nothing)

    # Compute drift at uprev: drift(X_n) = ν*a(X_n) into nlsolver.z (will be overwritten)
    c(nlsolver.z, uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up nlsolver
    # tmp = X_n + ν*k - θ*dt*drift(X_n)
    @.. nlsolver.tmp = uprev + nlsolver.tmp - theta * dt * nlsolver.z
    # Initial guess: z = dt * drift(uprev) (linear extrapolation)
    @.. nlsolver.z = dt * nlsolver.z

    # Step 4: Solve using nlsolver
    # The nlsolver solves: z = dt * drift(tmp + θ*z)
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    # Step 5: Final update: X_{n+1} = tmp + θ*z
    @.. u = nlsolver.tmp + theta * z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.rate(P.cache.currate, u, p, t + dt)
    end

    return nothing
end

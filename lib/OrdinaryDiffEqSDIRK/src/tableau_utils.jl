# Tableau-based utility functions for SDIRK methods

@inline function compute_sdirk_stage!(integrator, cache, nlsolver, stage::Int, 
                                     z_prev, stage_coeffs, c_val, tmp_val)
    nlsolver.z = z_prev
    nlsolver.tmp = tmp_val
    nlsolver.c = c_val
    z_new = nlsolve!(nlsolver, integrator, cache, false)
    nlsolvefail(nlsolver) && return nothing
    z_new
end

@inline function compute_stage_constantcache!(integrator, cache, stage::Int,
                                            prev_z, coeffs, c_val, base_tmp)
    nlsolver = cache.nlsolver
    nlsolver.z = prev_z
    nlsolver.tmp = base_tmp
    nlsolver.c = c_val
    z = nlsolve!(nlsolver, integrator, cache, false)
    nlsolvefail(nlsolver) && return nothing
    z
end

@inline function compute_stage_mutablecache!(integrator, cache, stage::Int,
                                           prev_z, coeffs, c_val, base_tmp)
    nlsolver = cache.nlsolver
    nlsolver.z = prev_z
    nlsolver.tmp = base_tmp  
    nlsolver.c = c_val
    z = nlsolve!(nlsolver, integrator, cache, false)
    nlsolvefail(nlsolver) && return nothing
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z
end

# generic error estimation for embedded methods
@inline function compute_embedded_error!(integrator, cache, btilde_coeffs, z_stages)
    if integrator.opts.adaptive
        tmp = sum(btilde_coeffs[i] * z_stages[i] for i in eachindex(z_stages))
        alg = unwrap_alg(integrator, true)
        nlsolver = cache.nlsolver
        
        if isnewton(nlsolver) && alg.smooth_est
            integrator.stats.nsolve += 1
            if hasfield(typeof(cache), :atmp)
                est = cache.atmp
                linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                    linu = _vec(est))
            else
                est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
            end
        else
            est = tmp
        end
        
        if hasfield(typeof(cache), :atmp)
            calculate_residuals!(cache.atmp, est, integrator.uprev, integrator.u, 
                               integrator.opts.abstol, integrator.opts.reltol, 
                               integrator.opts.internalnorm, integrator.t)
            integrator.EEst = integrator.opts.internalnorm(cache.atmp, integrator.t)
        else
            atmp = calculate_residuals(est, integrator.uprev, integrator.u, 
                                     integrator.opts.abstol, integrator.opts.reltol, 
                                     integrator.opts.internalnorm, integrator.t)
            integrator.EEst = integrator.opts.internalnorm(atmp, integrator.t)
        end
    end
end

# simple 1-stage SDIRK perform step (ImplicitEuler, ImplicitMidpoint)
@muladd function perform_step_1stage_sdirk!(integrator, cache, γ_val, adaptive_error_est=true, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else
        nlsolver.z = zero(u)
    end

    nlsolver.tmp = uprev
    nlsolver.γ = γ_val
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + z

    # error estimation for ImplicitEuler only
    if adaptive_error_est && integrator.opts.adaptive && integrator.success_iter > 0 && γ_val == 1.0
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
        
        if hasfield(typeof(cache), :atmp)
            tmp = nlsolver.tmp
            @.. broadcast=false tmp=r * integrator.opts.internalnorm(
                (u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
            calculate_residuals!(cache.atmp, tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.atmp, t)
        else
            tmp = r * integrator.opts.internalnorm.((u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
            atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    else
        integrator.EEst = 1
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    
    if hasfield(typeof(cache), :step_limiter!)
        cache.step_limiter!(u, integrator, p, t + dt)
    end
    
    if hasfield(typeof(cache), :atmp) && integrator.opts.adaptive && integrator.differential_vars !== nothing
        @.. broadcast=false cache.atmp=ifelse(cache.algebraic_vars, integrator.fsallast, false) /
                                 integrator.opts.abstol
        integrator.EEst += integrator.opts.internalnorm(cache.atmp, t)
    end
end

# mutable cache version of 1-stage SDIRK  
@muladd function perform_step_1stage_sdirk_mutable!(integrator, cache, γ_val, adaptive_error_est=true, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack nlsolver, step_limiter! = cache
    @unpack z, tmp = nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    if alg.extrapolant == :linear
        @.. broadcast=false z=dt * integrator.fsalfirst
    else
        z .= zero(eltype(u))
    end

    nlsolver.tmp .= uprev
    nlsolver.γ = γ_val
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=uprev + z

    step_limiter!(u, integrator, p, t + dt)

    # Error estimation for ImplicitEuler only
    if adaptive_error_est && integrator.opts.adaptive && integrator.success_iter > 0 && γ_val == 1.0
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
        @.. broadcast=false tmp=r * integrator.opts.internalnorm(
            (u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
        calculate_residuals!(cache.atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(cache.atmp, t)
    else
        integrator.EEst = 1
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(integrator.fsallast, u, p, t + dt)

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        @.. broadcast=false cache.atmp=ifelse(cache.algebraic_vars, integrator.fsallast, false) /
                                 integrator.opts.abstol
        integrator.EEst += integrator.opts.internalnorm(cache.atmp, t)
    end
end

# Generic 2-stage SDIRK for SDIRK2 and similar methods
@muladd function perform_step_2stage_sdirk!(integrator, cache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack γ, d, α1, α2, btilde1, btilde2 = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    z₁ = dt * integrator.fsalfirst
    nlsolver.z = z₁
    nlsolver.tmp = uprev + d * z₁
    nlsolver.c = γ
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    z₂ = α1 * z₁
    nlsolver.z = z₂
    nlsolver.tmp = uprev + d * z₁
    nlsolver.c = 1
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + d * z₂

    # Error estimation
    if integrator.opts.adaptive
        tmp = btilde1 * z₁ + btilde2 * z₂
        if isnewton(nlsolver) && alg.smooth_est
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₂ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

# Mutable cache version of 2-stage SDIRK
@muladd function perform_step_2stage_sdirk_mutable!(integrator, cache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, z₂, atmp, nlsolver, step_limiter! = cache
    @unpack tmp = nlsolver
    @unpack γ, d, α1, α2, btilde1, btilde2 = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    @.. broadcast=false z₁=dt * integrator.fsalfirst
    nlsolver.z = z₁
    @.. broadcast=false tmp=uprev + d * z₁
    nlsolver.c = γ
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    @.. broadcast=false z₂=α1 * z₁
    nlsolver.z = z₂
    @.. broadcast=false tmp=uprev + d * z₁
    nlsolver.c = 1
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + d * z₂
    step_limiter!(u, integrator, p, t + dt)

    # Error estimation
    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂
        if isnewton(nlsolver) && alg.smooth_est
            integrator.stats.nsolve += 1
            est = atmp
            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est))
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₂ / dt
end 
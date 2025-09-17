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




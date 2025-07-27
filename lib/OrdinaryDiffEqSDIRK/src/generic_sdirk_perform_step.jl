# Generic tableau-based perform_step! implementation for SDIRK methods

using OrdinaryDiffEqCore: @unpack, unwrap_alg, OrdinaryDiffEqCore, calculate_residuals
using OrdinaryDiffEqNonlinearSolve: markfirststage!, nlsolve!, nlsolvefail, isnewton, set_new_W!, get_W
using OrdinaryDiffEqDifferentiation: dolinsolve
using LinearAlgebra: I, mul!
using FastBroadcast: @..
using MuladdMacro: @muladd

function initialize!(integrator, cache::SDIRKConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::SDIRKMutableCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function generic_sdirk_perform_step!(integrator, cache::SDIRKConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tab, nlsolver = cache
    alg = unwrap_alg(integrator, true)
    
    s = size(tab.A, 1)
    c = tab.c
    A = tab.A
    b = tab.b
    b_embed = tab.b_embed
    γ = tab.γ
    
    z = Vector{typeof(u)}(undef, s)
    markfirststage!(nlsolver)
    
    for i in 1:s
        stage_sum = uprev
        for j in 1:i-1
            stage_sum += A[i,j] * z[j]
        end
        
        if i == 1 && alg.extrapolant == :linear
            z_guess = dt * integrator.fsalfirst
        else
            z_guess = zero(u)
        end
        
        nlsolver.z = z_guess
        nlsolver.tmp = stage_sum
        nlsolver.c = typeof(nlsolver.c)(c[i])
        nlsolver.γ = typeof(nlsolver.γ)(A[i,i])
        
        z[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
    end
    
    u = uprev
    for i in 1:s
        u += b[i] * z[i]
    end
    
    if integrator.opts.adaptive && b_embed !== nothing
        tmp = zero(u)
        for i in 1:s
            tmp += b_embed[i] * z[i]
        end
        
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
    
    integrator.fsallast = z[s] ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SDIRK2Cache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, z₂, atmp, nlsolver, tab, step_limiter! = cache
    @unpack tmp = nlsolver
    alg = unwrap_alg(integrator, true)
    
    A = tab.A
    b = tab.b
    c = tab.c
    b_embed = tab.b_embed
    
    markfirststage!(nlsolver)
    
    if alg.extrapolant == :linear
        @.. broadcast=false z₁ = dt * integrator.fsalfirst
    else
        fill!(z₁, zero(eltype(u)))
    end
    
    nlsolver.z = z₁
    @.. broadcast=false nlsolver.tmp = uprev
            nlsolver.c = typeof(nlsolver.c)(c[1])
        nlsolver.γ = typeof(nlsolver.γ)(A[1,1])
    
    z₁ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    
    @.. broadcast=false z₂ = zero(eltype(u))
    nlsolver.z = z₂
    @.. broadcast=false nlsolver.tmp = uprev + A[2,1] * z₁
            nlsolver.c = typeof(nlsolver.c)(c[2])
        nlsolver.γ = typeof(nlsolver.γ)(A[2,2])
    
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z₂ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    
    @.. broadcast=false u = uprev + b[1] * z₁ + b[2] * z₂
    step_limiter!(u, integrator, p, t + dt)
    
    if integrator.opts.adaptive && b_embed !== nothing
        @.. broadcast=false tmp = b_embed[1] * z₁ + b_embed[2] * z₂
        
        if alg.smooth_est && isnewton(nlsolver)
            est = atmp
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
    
    @.. broadcast=false integrator.fsallast = z₂ / dt
end

# Specialized dispatcher for each method that uses the generic implementation
# These just need to ensure the cache has the correct tableau

# Basic SDIRK methods
@muladd function perform_step!(integrator, cache::ImplicitEulerConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::TrapezoidConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SDIRK2ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SDIRK22ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SDIRK22Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::Cash4ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::Cash4Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK4ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK4Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK5ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK5Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK6ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK6Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK7ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK7Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK8ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SFSDIRK8Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK54I8L2SAConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK54I8L2SACache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK436L2SA2ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK436L2SA2Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK437L2SAConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK437L2SACache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK547L2SA2ConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK547L2SA2Cache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK659L2SAConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::ESDIRK659L2SACache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

# TRBDF2 special handling
@muladd function perform_step!(integrator, cache::TRBDF2ConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    tab = cache.tab
    γ = tab.c[2]
    d = tab.γ
    ω = tab.b[1]
    btilde1, btilde2, btilde3 = tab.b_embed[1], tab.b_embed[2], tab.b_embed[3]
    α1, α2 = tab.A[3,1], tab.A[3,2]
    
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    zprev = dt * integrator.fsalfirst

    zγ = zprev
    nlsolver.z = zγ
    nlsolver.c = γ
    nlsolver.tmp = uprev + d * zprev
    zγ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    z = α1 * zprev + α2 * zγ
    nlsolver.z = z
    nlsolver.c = 1
    nlsolver.tmp = uprev + ω * zprev + ω * zγ
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + d * z

    if integrator.opts.adaptive
        tmp = btilde1 * zprev + btilde2 * zγ + btilde3 * z
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

    integrator.fsallast = z ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack zprev, zγ, atmp, nlsolver, step_limiter! = cache
    @unpack z, tmp = nlsolver
    tab = cache.tab
    γ = tab.c[2]
    d = tab.γ
    ω = tab.b[1]
    btilde1, btilde2, btilde3 = tab.b_embed[1], tab.b_embed[2], tab.b_embed[3]
    α1, α2 = tab.A[3,1], tab.A[3,2]
    
    alg = unwrap_alg(integrator, true)

    @.. broadcast=false zprev = dt * integrator.fsalfirst
    markfirststage!(nlsolver)

    @.. broadcast=false zγ = zprev
    z .= zγ
    @.. broadcast=false tmp = uprev + d * zprev
    nlsolver.c = γ
    zγ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false z = α1 * zprev + α2 * zγ
    @.. broadcast=false tmp = uprev + ω * zprev + ω * zγ
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u = tmp + d * z
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false tmp = btilde1 * zprev + btilde2 * zγ + btilde3 * z
        if alg.smooth_est && isnewton(nlsolver)
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

    @.. broadcast=false integrator.fsallast = z / dt
end



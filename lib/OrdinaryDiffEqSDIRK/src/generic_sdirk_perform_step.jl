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
    # apply step limiter if available on algorithm
    if hasproperty(alg, :step_limiter!)
        alg.step_limiter!(u, integrator, p, t + dt)
    end
    
    if integrator.opts.adaptive && b_embed !== nothing
        tmp = zero(u)
        for i in 1:s
            tmp += b_embed[i] * z[i]
        end
        
        has_smooth_est = hasfield(typeof(alg), :smooth_est)
        if isnewton(nlsolver) && has_smooth_est && alg.smooth_est
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

@muladd function generic_sdirk_perform_step!(integrator, cache::SDIRKMutableCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack z₁, atmp, nlsolver, tab, step_limiter! = cache
    @unpack tmp = nlsolver
    alg = unwrap_alg(integrator, true)
    
    A = tab.A
    b = tab.b
    c = tab.c
    b_embed = tab.b_embed
    s = size(A, 1)
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
    
    # handle additional stages if present
    z₂ = nothing
    if s >= 2
        if hasfield(typeof(cache), :z₂)
            z₂ = getfield(cache, :z₂)
            @.. broadcast=false z₂ = zero(eltype(u))
            nlsolver.z = z₂
            @.. broadcast=false nlsolver.tmp = uprev + A[2,1] * z₁
            nlsolver.c = typeof(nlsolver.c)(c[2])
            nlsolver.γ = typeof(nlsolver.γ)(A[2,2])
            
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
            z₂ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
        end
    end
    
    # computes final solution
    if s == 1
        @.. broadcast=false u = uprev + b[1] * z₁
    elseif s == 2 && z₂ !== nothing
        @.. broadcast=false u = uprev + b[1] * z₁ + b[2] * z₂
    else
        error("Unsupported number of stages: $s")
    end
    
    step_limiter!(u, integrator, p, t + dt)
    
    # error estimation for adaptive methods
    if integrator.opts.adaptive && b_embed !== nothing
        if s == 1
            @.. broadcast=false tmp = b_embed[1] * z₁
        elseif s == 2 && z₂ !== nothing
            @.. broadcast=false tmp = b_embed[1] * z₁ + b_embed[2] * z₂
        end
        
        has_smooth_est = hasfield(typeof(alg), :smooth_est)
        if has_smooth_est && alg.smooth_est && isnewton(nlsolver)
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
    
    # set fsallast based on last stage and update integrator state
    if s == 1
        @.. broadcast=false integrator.fsallast = z₁ / dt
    elseif s == 2 && z₂ !== nothing
        @.. broadcast=false integrator.fsallast = z₂ / dt
    end
    # keep k array consistent with FSAL bookkeeping
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function generic_additive_sdirk_perform_step!(integrator, cache::SDIRKConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tab, nlsolver = cache
    alg = unwrap_alg(integrator, true)
    
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f_expl = integrator.f.f2
    else
        f_impl = integrator.f
    end
    
    markfirststage!(nlsolver)
    
    s = size(tab.A, 1)
    γ = tab.γ
    A = tab.A
    b = tab.b
    c = tab.c
    b_embed = tab.b_embed
    
    A_explicit = tab.A_explicit
    b_explicit = tab.b_explicit
    
    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        z₁ = dt * f_impl(uprev, p, t)
    else
        z₁ = dt * integrator.fsalfirst
    end
    
    if s >= 2 && hasfield(typeof(cache), :z₂)
        if hasproperty(tab, :α_pred) && tab.α_pred !== nothing
            z₂ = getfield(cache, :z₂)
            @.. broadcast=false z₂ = zero(eltype(u))
            for j in 1:1
                @.. broadcast=false z₂ += tab.α_pred[2, j] * (j == 1 ? z₁ : zero(z₁))
            end
        else
            # Fallback: copy previous stage
            z₂ = z₁
        end
        nlsolver.z = z₂
        tmp = uprev + γ * z₁
        
        if integrator.f isa SplitFunction
            k1 = dt * integrator.fsalfirst - z₁
            tmp += A_explicit[2,1] * k1
        end
        
        nlsolver.tmp = tmp
        nlsolver.c = 2γ
        z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        
        if integrator.f isa SplitFunction
            z₃ = z₂
            u = nlsolver.tmp + γ * z₂
            k2 = dt * f_expl(u, p, t + 2γ * dt)
            integrator.stats.nf2 += 1
            tmp = uprev + A[3,1] * z₁ + A[3,2] * z₂ + A_explicit[3,1] * k1 + A_explicit[3,2] * k2
        else
            if hasproperty(tab, :α_pred) && tab.α_pred !== nothing
                z₃ = tab.α_pred[3,1] * z₁ + tab.α_pred[3,2] * z₂
            else
                θ = c[3] / c[2]
                α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c[2]) * γ)
                α32 = ((-2θ + 3θ^2)) + (6θ * (1 - θ) / c[2]) * γ
                z₃ = α31 * z₁ + α32 * z₂
            end
            tmp = uprev + A[3,1] * z₁ + A[3,2] * z₂
        end
        
        nlsolver.z = z₃
        nlsolver.tmp = tmp
        nlsolver.c = c[3]
        z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        
        # Step 4
        if integrator.f isa SplitFunction
            z₄ = z₂
            u = nlsolver.tmp + γ * z₃
            k3 = dt * f_expl(u, p, t + c[3] * dt)
            integrator.stats.nf2 += 1
            tmp = uprev + A[4,1] * z₁ + A[4,2] * z₂ + A[4,3] * z₃ + A_explicit[4,1] * k1 + A_explicit[4,2] * k2 + A_explicit[4,3] * k3
        else
            if hasproperty(tab, :α_pred) && tab.α_pred !== nothing
                z₄ = tab.α_pred[4,1] * z₁ + tab.α_pred[4,2] * z₂ + tab.α_pred[4,3] * z₃
            else
                z₄ = z₁
            end
            tmp = uprev + A[4,1] * z₁ + A[4,2] * z₂ + A[4,3] * z₃
        end
        
        nlsolver.z = z₄
        nlsolver.tmp = tmp
        nlsolver.c = 1
        z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        
        u = nlsolver.tmp + γ * z₄
        if integrator.f isa SplitFunction
            k4 = dt * f_expl(u, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + A[4,1] * z₁ + A[4,2] * z₂ + A[4,3] * z₃ + γ * z₄ + 
                b_explicit[1] * k1 + b_explicit[2] * k2 + b_explicit[3] * k3 + b_explicit[4] * k4
        end
        
        # Error estimation
        if integrator.opts.adaptive && b_embed !== nothing
            if integrator.f isa SplitFunction
                tmp = b_embed[1] * z₁ + b_embed[2] * z₂ + b_embed[3] * z₃ + b_embed[4] * z₄
            else
                tmp = b_embed[1] * z₁ + b_embed[2] * z₂ + b_embed[3] * z₃ + b_embed[4] * z₄
            end
            
            has_smooth_est = hasfield(typeof(alg), :smooth_est)
            if isnewton(nlsolver) && has_smooth_est && alg.smooth_est
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
        
    else
        generic_sdirk_perform_step!(integrator, cache, repeat_step)
    end
end

@muladd function perform_step!(integrator, cache::Union{
    SDIRK2ConstantCache, SDIRK2Cache,
    SFSDIRK4ConstantCache, SFSDIRK4Cache,
    SFSDIRK5ConstantCache, SFSDIRK5Cache,
    SFSDIRK6ConstantCache, SFSDIRK6Cache,
    SFSDIRK7ConstantCache, SFSDIRK7Cache,
    SFSDIRK8ConstantCache, SFSDIRK8Cache,
    ESDIRK54I8L2SAConstantCache, ESDIRK54I8L2SACache,
    ESDIRK436L2SA2ConstantCache, ESDIRK436L2SA2Cache,
    ESDIRK437L2SAConstantCache, ESDIRK437L2SACache,
    ESDIRK547L2SA2ConstantCache, ESDIRK547L2SA2Cache,
    ESDIRK659L2SAConstantCache, ESDIRK659L2SACache
}, repeat_step=false)
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
    alg.step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        tmp = btilde1 * zprev + btilde2 * zγ + btilde3 * z
        has_smooth_est = hasfield(typeof(alg), :smooth_est)
        if isnewton(nlsolver) && has_smooth_est && alg.smooth_est
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
    @unpack zprev, zᵧ, atmp, nlsolver, step_limiter! = cache
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

    @.. broadcast=false zᵧ = zprev
    z .= zᵧ
    @.. broadcast=false tmp = uprev + d * zprev
    nlsolver.c = γ
    zᵧ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false z = α1 * zprev + α2 * zᵧ
    @.. broadcast=false tmp = uprev + ω * zprev + ω * zᵧ
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u = tmp + d * z
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false tmp = btilde1 * zprev + btilde2 * zᵧ + btilde3 * z
        has_smooth_est = hasfield(typeof(alg), :smooth_est)
        if has_smooth_est && alg.smooth_est && isnewton(nlsolver)
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

# Hairer4 method
@muladd function perform_step!(integrator, cache::Union{Hairer4ConstantCache, Hairer4Cache}, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::Union{ImplicitEulerConstantCache, ImplicitEulerCache}, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack nlsolver, tab = cache
    alg = unwrap_alg(integrator, true)
    
    markfirststage!(nlsolver)
    
    γ = tab.A[1,1]
    
    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else
        nlsolver.z = zero(uprev)
    end
    
    nlsolver.tmp = uprev
    nlsolver.c = γ
    nlsolver.γ = γ
    
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    
    integrator.u = uprev + z
    # step limiter
    alg.step_limiter!(integrator.u, integrator, p, t + dt)
    integrator.fsallast = z ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack nlsolver, tab = cache
    alg = unwrap_alg(integrator, true)

    markfirststage!(nlsolver)

    γ = tab.A[1,1]
    c = tab.c[1]
    b1 = tab.b[1]

    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else
        nlsolver.z = zero(uprev)
    end

    nlsolver.tmp = uprev
    nlsolver.c = c
    nlsolver.γ = γ

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = uprev + b1 * z
    if hasproperty(alg, :step_limiter!)
        alg.step_limiter!(u, integrator, p, t + dt)
    end

    integrator.fsallast = z ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointCache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack z₁, nlsolver, tab, step_limiter! = cache
    @unpack tmp = nlsolver
    alg = unwrap_alg(integrator, true)

    γ = tab.A[1,1]
    c = tab.c[1]
    b1 = tab.b[1]

    markfirststage!(nlsolver)

    if alg.extrapolant == :linear
        @.. broadcast=false z₁ = dt * integrator.fsalfirst
    else
        fill!(z₁, zero(eltype(u)))
    end

    nlsolver.z = z₁
    @.. broadcast=false nlsolver.tmp = uprev
    nlsolver.c = typeof(nlsolver.c)(c)
    nlsolver.γ = typeof(nlsolver.γ)(γ)

    z₁ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u = uprev + b1 * z₁
    step_limiter!(u, integrator, p, t + dt)

    @.. broadcast=false integrator.fsallast = z₁ / dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Union{TrapezoidConstantCache, TrapezoidCache}, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::Union{SDIRK22ConstantCache, SDIRK22Cache}, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::Union{SSPSDIRK2ConstantCache, SSPSDIRK2Cache}, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end


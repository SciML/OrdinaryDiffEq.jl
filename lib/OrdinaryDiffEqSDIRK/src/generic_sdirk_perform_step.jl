# Generic tableau-based perform_step! implementation for SDIRK methods

using OrdinaryDiffEqCore: @unpack, unwrap_alg, OrdinaryDiffEqCore, calculate_residuals
using OrdinaryDiffEqNonlinearSolve: markfirststage!, nlsolve!, nlsolvefail, isnewton, set_new_W!, get_W
using OrdinaryDiffEqDifferentiation: dolinsolve
using LinearAlgebra: I, mul!
using FastBroadcast: @..
using MuladdMacro: @muladd

@inline _get_step_limiter(alg, cache) = hasproperty(alg, :step_limiter!) ? alg.step_limiter! : trivial_limiter!

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
    
    # Check if this is an IMEX (split function) scheme
    is_imex = integrator.f isa SplitFunction && tab.has_additive_splitting
    
    z = Vector{typeof(u)}(undef, s)
    
    # For IMEX schemes, we need to store explicit stage derivatives
    if is_imex
        k_explicit = Vector{typeof(u)}(undef, s)
        f_impl = integrator.f.f1
        f_expl = integrator.f.f2
    else
        k_explicit = nothing
        f_impl = integrator.f
        f_expl = nothing
    end
    
    markfirststage!(nlsolver)
    
    for i in 1:s
        # Compute implicit stage sum
        stage_sum = uprev
        for j in 1:i-1
            stage_sum += A[i,j] * z[j]
        end
        
        # For IMEX schemes, add explicit contributions
        if is_imex && tab.A_explicit !== nothing
            # First stage for explicit part
            if i == 1
                k_explicit[1] = dt * integrator.fsalfirst - (is_imex ? dt * f_impl(uprev, p, t) : zero(u))
            else
                # Compute intermediate solution for explicit evaluation
                u_tmp = nlsolver.tmp
                if i >= 2
                    u_tmp += tab.γ * z[i-1]
                end
                
                # Use c from implicit tableau if c_explicit is not defined
                c_exp = tab.c_explicit !== nothing ? tab.c_explicit[i] : tab.c[i]
                k_explicit[i] = dt * f_expl(u_tmp, p, t + c_exp * dt)
                integrator.stats.nf2 += 1
            end
            
            # Add explicit tableau contributions
            for j in 1:i-1
                stage_sum += tab.A_explicit[i,j] * k_explicit[j]
            end
        end

        # Determine initial guess for stage
        if i == 1
            if is_imex
                z_guess = dt * f_impl(uprev, p, t)
            else
                z_guess = (alg.extrapolant == :linear) ? dt * integrator.fsalfirst : zero(u)
            end
        elseif i > 1 && hasproperty(tab, :α_pred) && tab.α_pred !== nothing
            z_guess = zero(u)
            @inbounds for j in 1:i-1
                z_guess += tab.α_pred[i,j] * z[j]
            end
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
    
    # Compute final solution
    u = uprev
    for i in 1:s
        u += b[i] * z[i]
    end
    
    # For IMEX schemes, add explicit contributions to final solution
    if is_imex && tab.b_explicit !== nothing
        # Need final explicit evaluation
        if s >= 1
            u_final = nlsolver.tmp + tab.γ * z[s]
            k_explicit[s] = dt * f_expl(u_final, p, t + dt)
            integrator.stats.nf2 += 1
        end
        
        for i in 1:s
            u += tab.b_explicit[i] * k_explicit[i]
        end
    end
    
    # apply step limiter if available on algorithm
    if hasproperty(alg, :step_limiter!)
        alg.step_limiter!(u, integrator, p, t + dt)
    end
    
    # Error estimation
    if integrator.opts.adaptive && b_embed !== nothing
        tmp = zero(u)
        for i in 1:s
            tmp += b_embed[i] * z[i]
        end
        
        # For IMEX schemes, include explicit error contributions
        # (This would require additional error estimation tableau components)
        
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
    
    # Set final derivative appropriately
    if is_imex
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z[s] ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function generic_sdirk_perform_step!(integrator, cache::SDIRKMutableCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack zs, atmp, nlsolver, tab = cache
    @unpack tmp = nlsolver
    alg = unwrap_alg(integrator, true)
    step_limiter! = _get_step_limiter(alg, cache)
    
    A = tab.A
    b = tab.b
    c = tab.c
    b_embed = tab.b_embed
    s = size(A, 1)
    
    # Check if this is an IMEX (split function) scheme
    is_imex = integrator.f isa SplitFunction && tab.has_additive_splitting
    
    # For IMEX schemes, we need to store explicit stage derivatives
    if is_imex
        k_explicit = Vector{typeof(u)}(undef, s)
        for i in 1:s
            k_explicit[i] = zero(u)
        end
        f_impl = integrator.f.f1
        f_expl = integrator.f.f2
    else
        k_explicit = nothing
        f_impl = integrator.f
        f_expl = nothing
    end
    
    markfirststage!(nlsolver)
    
    # Solve each stage
    for i in 1:s
        zi = zs[i]
        
        # Determine initial guess for stage
        if i == 1
            if is_imex && !repeat_step && !integrator.last_stepfail
                # For IMEX, explicit tableau is not FSAL
                f_impl(zi, uprev, p, t)
                @.. broadcast=false zi *= dt
            elseif alg.extrapolant == :linear
                @.. broadcast=false zi = dt * integrator.fsalfirst
            else
                fill!(zi, zero(eltype(u)))
            end
        else
            # Initialize stage to zero
            fill!(zi, zero(eltype(u)))
            
            # Add predictor if available
            if hasproperty(tab, :α_pred) && tab.α_pred !== nothing
                for j in 1:i-1
                    @.. broadcast=false zi += tab.α_pred[i, j] * zs[j]
                end
            end
        end
        
        nlsolver.z = zi
        
        # Compute implicit stage sum
        @.. broadcast=false nlsolver.tmp = uprev
        for j in 1:i-1
            @.. broadcast=false nlsolver.tmp += A[i, j] * zs[j]
        end
        
        # For IMEX schemes, add explicit contributions
        if is_imex && tab.A_explicit !== nothing
            # First stage for explicit part
            if i == 1
                @.. broadcast=false k_explicit[1] = dt * integrator.fsalfirst - zi
            else
                # Compute intermediate solution for explicit evaluation
                @.. broadcast=false u = nlsolver.tmp + A[i,i] * zs[i-1]
                # Use c from implicit tableau if c_explicit is not defined
                c_exp = tab.c_explicit !== nothing ? tab.c_explicit[i] : tab.c[i]
                f_expl(k_explicit[i], u, p, t + c_exp * dt)
                @.. broadcast=false k_explicit[i] *= dt
                integrator.stats.nf2 += 1
            end
            
            # Add explicit tableau contributions
            for j in 1:i-1
                @.. broadcast=false nlsolver.tmp += tab.A_explicit[i, j] * k_explicit[j]
            end
        end
        
        nlsolver.c = typeof(nlsolver.c)(c[i])
        nlsolver.γ = typeof(nlsolver.γ)(A[i, i])
        
        if i > 1 && isnewton(nlsolver)
            set_new_W!(nlsolver, false)
        end
        
        zi .= nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
    end
    
    # Compute final solution
    @.. broadcast=false u = uprev
    for i in 1:s
        @.. broadcast=false u += b[i] * zs[i]
    end
    
    # For IMEX schemes, add explicit contributions to final solution
    if is_imex && tab.b_explicit !== nothing
        # Need final explicit evaluation
        if s >= 1
            @.. broadcast=false tmp = nlsolver.tmp + A[s,s] * zs[s]
            f_expl(k_explicit[s], tmp, p, t + dt)
            @.. broadcast=false k_explicit[s] *= dt
            integrator.stats.nf2 += 1
        end
        
        for i in 1:s
            @.. broadcast=false u += tab.b_explicit[i] * k_explicit[i]
        end
    end
    
    step_limiter!(u, integrator, p, t + dt)
    
    # Error estimation for adaptive methods
    if integrator.opts.adaptive && b_embed !== nothing
        @.. broadcast=false tmp = zero(eltype(u))
        for i in 1:s
            @.. broadcast=false tmp += b_embed[i] * zs[i]
        end
        
        # For IMEX schemes, include explicit error contributions
        # (This would require additional error estimation tableau components)
        
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
    
    # Set final derivative appropriately
    if is_imex
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast = zs[s] / dt
    end
    
    # Keep k array consistent with FSAL bookkeeping
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

# Dispatch for unified caches - all SDIRK algorithms use these same cache types
@muladd function perform_step!(integrator, cache::SDIRKCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
end

@muladd function perform_step!(integrator, cache::SDIRKConstantCache, repeat_step=false)
    generic_sdirk_perform_step!(integrator, cache, repeat_step)
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

# TRBDF2 special handling for constant cache
@muladd function perform_step_trbdf2!(integrator, cache::SDIRKConstantCache, repeat_step=false)
    @unpack t, dt, uprev, u, f, p = integrator
    tab = cache.tab
    γ = tab.c[2]
    d = tab.γ
    ω = tab.b[1]
    btilde1, btilde2, btilde3 = tab.b_embed[1], tab.b_embed[2], tab.b_embed[3]
    α1, α2 = tab.A[3,1], tab.A[3,2]
    
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    step_limiter! = _get_step_limiter(alg, cache)
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
    step_limiter!(u, integrator, p, t + dt)

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

    # TRBDF2 is stiffly accurate; use fresh f at t+dt for FSAL
    integrator.fsallast = integrator.f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

# TRBDF2 special handling for mutable cache 
@muladd function perform_step_trbdf2!(integrator, cache::SDIRKCache, repeat_step=false)
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

    integrator.fsallast = integrator.f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end



# All other SDIRK methods use the generic implementation through the unified cache dispatch above


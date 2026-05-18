function initialize!(integrator, cache::ESDIRKIMEXConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function initialize!(integrator, cache::ESDIRKIMEXCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ESDIRKIMEXConstantCache, repeat_step = false
    )
    (; t, dt, uprev, u, p) = integrator
    nlsolver = cache.nlsolver
    tab = cache.tab
    (; Ai, bi, Ae, be, c, btilde, ebtilde, α, s) = tab
    alg = unwrap_alg(integrator, true)
    γ = Ai[s, s]

    f2 = nothing
    k = Vector{typeof(u)}(undef, s)
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    z = Vector{typeof(u)}(undef, s)

    markfirststage!(nlsolver)

    if tab.explicit_first_stage
        if integrator.f isa SplitFunction
            z[1] = dt * f_impl(uprev, p, t)
        else
            z[1] = dt * integrator.fsalfirst
        end
        if integrator.f isa SplitFunction
            k[1] = dt * integrator.fsalfirst - z[1]
        end
    else
        # Implicit first stage: Ai[1,1] ≠ 0, requires an nlsolve
        if integrator.success_iter > 0 && !integrator.reeval_fsal &&
                alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                alg.extrapolant == :interpolant
            current_extrapolant!(u, t + dt, integrator)
            z[1] = u - uprev
        elseif tab.stage1_extrapolation &&
                alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                alg.extrapolant == :linear
            z[1] = dt * integrator.fsalfirst
        else
            z[1] = zero(u)
        end
        nlsolver.z = z[1]
        nlsolver.tmp = uprev
        nlsolver.c = c[1]
        nlsolver.γ = γ
        z[1] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
    end

    for i in 2:s
        tmp = uprev
        for j in 1:(i - 1)
            tmp = tmp + Ai[i, j] * z[j]
        end

        if integrator.f isa SplitFunction
            for j in 1:(i - 1)
                tmp = tmp + Ae[i, j] * k[j]
            end
        end

        if integrator.f isa SplitFunction
            z_guess = z[1]
        elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[i])
            z_guess = tab.const_stage_guess[i]
        elseif !isempty(α) && !iszero(α[i])
            z_guess = zero(u)
            for j in 1:(i - 1)
                z_guess = z_guess + α[i][j] * z[j]
            end
        else
            z_guess = zero(u)
        end

        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[i]
        nlsolver.γ = γ
        z[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return

        if integrator.f isa SplitFunction && i < s
            u_stage = tmp + γ * z[i]
            k[i] = dt * f2(u_stage, p, t + c[i] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if integrator.f isa SplitFunction
        u = nlsolver.tmp + γ * z[s]
        k[s] = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev
        for i in 1:s
            u = u + bi[i] * z[i] + be[i] * k[i]
        end
    elseif tab.stiffly_accurate
        # b == A[s,:] (stiffly accurate) ⟹ u = uprev + Σbᵢzᵢ = tmp + γ*z[s] (Hairer & Wanner II, §IV.8)
        u = nlsolver.tmp + γ * z[s]
    else
        u = uprev
        for i in 1:s
            u = u + bi[i] * z[i]
        end
    end

    if integrator.opts.adaptive && !isempty(btilde)
        tmp = zero(u)
        for i in 1:s
            tmp = tmp + btilde[i] * z[i]
        end
        if integrator.f isa SplitFunction && !isempty(ebtilde)
            for i in 1:s
                tmp = tmp + ebtilde[i] * k[i]
            end
        end
        if isnewton(nlsolver) && _esdirk_smooth_est(alg)
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(
            est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    elseif tab.explicit_fsallast
        integrator.fsallast = integrator.f(u, p, t + tab.fsallast_c * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z[s] ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ESDIRKIMEXCache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; zs, ks, atmp, nlsolver, step_limiter!) = cache
    (; tmp) = nlsolver
    tab = cache.tab
    (; Ai, bi, Ae, be, c, btilde, ebtilde, α, s, reuse_W_at_stage2, split_guess) = tab
    alg = unwrap_alg(integrator, true)
    γ = Ai[s, s]

    f2 = nothing
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    if tab.explicit_first_stage
        if integrator.f isa SplitFunction && tab.fsal && !repeat_step && !integrator.last_stepfail
            f_impl(zs[1], integrator.uprev, p, integrator.t)
            zs[1] .*= dt
        else
            @..zs[1] = dt * integrator.fsalfirst
        end

        if integrator.f isa SplitFunction
            @..ks[1] = dt * integrator.fsalfirst - zs[1]
        end

        for i in 2:s
            copyto!(tmp, uprev)
            for j in 1:(i - 1)
                @..tmp = tmp + Ai[i, j] * zs[j]
            end

            if integrator.f isa SplitFunction
                for j in 1:(i - 1)
                    @..tmp = tmp + Ae[i, j] * ks[j]
                end
            end

            if integrator.f isa SplitFunction && split_guess[i] > 0
                copyto!(zs[i], zs[split_guess[i]])
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[i])
                fill!(zs[i], tab.const_stage_guess[i])
            elseif !isempty(α) && !iszero(α[i])
                fill!(zs[i], zero(eltype(u)))
                for j in 1:(i - 1)
                    @..zs[i] = zs[i] + α[i][j] * zs[j]
                end
            else
                fill!(zs[i], zero(eltype(u)))
            end

            nlsolver.z = zs[i]
            nlsolver.tmp = tmp
            nlsolver.c = c[i]
            zs[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            if i == 2 && reuse_W_at_stage2
                isnewton(nlsolver) && set_new_W!(nlsolver, false)
            end

            if integrator.f isa SplitFunction && i < s
                @..u = tmp + γ * zs[i]
                f2(ks[i], u, p, t + c[i] * dt)
                ks[i] .*= dt
                integrator.stats.nf2 += 1
            end
        end
    else
        # Implicit first stage: Ai[1,1] ≠ 0, requires an nlsolve
        if integrator.success_iter > 0 && !integrator.reeval_fsal &&
                alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                alg.extrapolant == :interpolant
            current_extrapolant!(u, t + dt, integrator)
            @.. broadcast = false zs[1] = u - uprev
        elseif tab.stage1_extrapolation &&
                alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                alg.extrapolant == :linear
            @.. broadcast = false zs[1] = dt * integrator.fsalfirst
        else
            zs[1] .= zero(eltype(zs[1]))
        end
        nlsolver.z = zs[1]
        nlsolver.tmp = uprev
        zs[1] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        # All stages share the same W (constant diagonal γ); reuse from stage 1
        isnewton(nlsolver) && set_new_W!(nlsolver, false)

        for i in 2:s
            copyto!(tmp, uprev)
            for j in 1:(i - 1)
                @..tmp = tmp + Ai[i, j] * zs[j]
            end

            if !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[i])
                fill!(zs[i], tab.const_stage_guess[i])
            elseif !isempty(α) && !iszero(α[i])
                fill!(zs[i], zero(eltype(u)))
                for j in 1:(i - 1)
                    @..zs[i] = zs[i] + α[i][j] * zs[j]
                end
            else
                fill!(zs[i], zero(eltype(u)))
            end

            nlsolver.z = zs[i]
            nlsolver.tmp = tmp
            nlsolver.c = c[i]
            zs[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
        end
    end

    if integrator.f isa SplitFunction
        @..u = tmp + γ * zs[s]
        f2(ks[s], u, p, t + dt)
        ks[s] .*= dt
        integrator.stats.nf2 += 1
        copyto!(u, uprev)
        for i in 1:s
            @..u = u + bi[i] * zs[i] + be[i] * ks[i]
        end
    elseif tab.stiffly_accurate
        # b == A[s,:] (stiffly accurate) ⟹ u = uprev + Σbᵢzᵢ = tmp + γ*z[s] (Hairer & Wanner II, §IV.8)
        @..u = tmp + γ * zs[s]
    else
        copyto!(u, uprev)
        for i in 1:s
            @..u = u + bi[i] * zs[i]
        end
    end

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive && !isempty(btilde)
        fill!(tmp, zero(eltype(u)))
        for i in 1:s
            @..tmp = tmp + btilde[i] * zs[i]
        end
        if integrator.f isa SplitFunction && !isempty(ebtilde)
            for i in 1:s
                @..tmp = tmp + ebtilde[i] * ks[i]
            end
        end
        if isnewton(nlsolver) && _esdirk_smooth_est(alg)
            est = nlsolver.cache.dz
            linres = dolinsolve(
                integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est)
            )
            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(
            atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    elseif tab.explicit_fsallast
        integrator.f(integrator.fsallast, u, p, t + tab.fsallast_c * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    else
        @..integrator.fsallast = zs[s] / dt
    end
end

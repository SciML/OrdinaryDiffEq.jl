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
    _perform_step_oop!(integrator, cache, repeat_step, cache.tab)
end

@muladd function perform_step!(integrator, cache::ESDIRKIMEXCache, repeat_step = false)
    _perform_step_iip!(integrator, cache, repeat_step, cache.tab)
end

function _build_fused_accum(dst::Symbol, base, A::Symbol, zs::Symbol, i::Int, jrange)
    if isempty(jrange)
        return :(copyto!($dst, $base))
    end
    terms = [:($A[$i, $j] * $zs[$j]) for j in jrange]
    rhs = base === nothing ? Expr(:call, :+, terms...) : Expr(:call, :+, base, terms...)
    return :(@.. broadcast = false $dst = $rhs)
end

function _build_alpha_guess(i::Int)
    if i == 1
        return :(fill!(zs[$i], zero(eltype(u))))
    end
    terms = [:(α[$i][$j] * zs[$j]) for j in 1:(i - 1)]
    rhs = length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
    return :(@.. broadcast = false zs[$i] = $rhs)
end

function _build_initial_guess_efs_false(i::Int)
    αfused = _build_alpha_guess(i)
    return quote
        if !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
            fill!(zs[$i], tab.const_stage_guess[$i])
        elseif !isempty(α) && !iszero(α[$i])
            $αfused
        else
            fill!(zs[$i], zero(eltype(u)))
        end
    end
end

function _build_initial_guess_efs_true(i::Int)
    αfused = _build_alpha_guess(i)
    return quote
        if integrator.f isa SplitFunction && split_guess[$i] > 0
            copyto!(zs[$i], zs[split_guess[$i]])
        elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
            fill!(zs[$i], tab.const_stage_guess[$i])
        elseif !isempty(α) && !iszero(α[$i])
            $αfused
        else
            fill!(zs[$i], zero(eltype(u)))
        end
    end
end

function _build_stages_efs_true(S::Int)
    body = quote end
    for i in 2:S
        tmp_accum = _build_fused_accum(:tmp, :uprev, :Ai, :zs, i, 1:(i - 1))

        ks_accum = quote end
        for j in 1:(i - 1)
            push!(ks_accum.args, :(@.. broadcast = false tmp = tmp + Ae[$i, $j] * ks[$j]))
        end

        initial_guess = _build_initial_guess_efs_true(i)

        stage = quote
            $tmp_accum
            if integrator.f isa SplitFunction
                $ks_accum
            end
            $initial_guess
            nlsolver.z = zs[$i]
            nlsolver.tmp = tmp
            nlsolver.c = c[$i]
            zs[$i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            $(i == 2 ? :(if reuse_W_at_stage2; isnewton(nlsolver) && set_new_W!(nlsolver, false); end) : nothing)
            $(i < S ? quote
                if integrator.f isa SplitFunction
                    @.. broadcast = false u = tmp + γ * zs[$i]
                    f2(ks[$i], u, p, t + c[$i] * dt)
                    ks[$i] .*= dt
                    integrator.stats.nf2 += 1
                end
            end : nothing)
        end
        append!(body.args, stage.args)
    end
    return body
end

function _build_stages_efs_false(S::Int)
    body = quote end
    for i in 2:S
        tmp_accum = _build_fused_accum(:tmp, :uprev, :Ai, :zs, i, 1:(i - 1))
        initial_guess = _build_initial_guess_efs_false(i)
        stage = quote
            $tmp_accum
            $initial_guess
            nlsolver.z = zs[$i]
            nlsolver.tmp = tmp
            nlsolver.c = c[$i]
            zs[$i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
        end
        append!(body.args, stage.args)
    end
    return body
end

function _build_output(S::Int)
    u_terms = [:(bi[$i] * zs[$i]) for i in 1:S]
    u_rhs = Expr(:call, :+, :uprev, u_terms...)
    u_split_terms = vcat(u_terms, [:(be[$i] * ks[$i]) for i in 1:S])
    u_split_rhs = Expr(:call, :+, :uprev, u_split_terms...)
    return quote
        if integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[$S]
            f2(ks[$S], u, p, t + dt)
            ks[$S] .*= dt
            integrator.stats.nf2 += 1
            @.. broadcast = false u = $u_split_rhs
        elseif tab.stiffly_accurate
            @.. broadcast = false u = tmp + γ * zs[$S]
        else
            @.. broadcast = false u = $u_rhs
        end
    end
end

function _build_adaptive(S::Int)
    btilde_terms = [:(btilde[$i] * zs[$i]) for i in 1:S]
    btilde_rhs = Expr(:call, :+, btilde_terms...)
    ebtilde_terms = [:(ebtilde[$i] * ks[$i]) for i in 1:S]
    ebtilde_rhs = Expr(:call, :+, ebtilde_terms...)
    return quote
        if integrator.opts.adaptive && !isempty(btilde)
            @.. broadcast = false tmp = $btilde_rhs
            if integrator.f isa SplitFunction && !isempty(ebtilde)
                @.. broadcast = false tmp = tmp + $ebtilde_rhs
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
    end
end

@generated function _perform_step_iip!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{S, T, T2}
    ) where {S, T, T2}
    setup = quote
        (; t, dt, uprev, u, p) = integrator
        (; zs, ks, atmp, nlsolver, step_limiter!) = cache
        (; tmp) = nlsolver
        Ai = tab.Ai
        bi = tab.bi
        Ae = tab.Ae
        be = tab.be
        c = tab.c
        btilde = tab.btilde
        ebtilde = tab.ebtilde
        α = tab.α
        reuse_W_at_stage2 = tab.reuse_W_at_stage2
        split_guess = tab.split_guess
        alg = unwrap_alg(integrator, true)
        γ = Ai[$S, $S]

        f2 = nothing
        if integrator.f isa SplitFunction
            f_impl = integrator.f.f1
            f2 = integrator.f.f2
        else
            f_impl = integrator.f
        end

        markfirststage!(nlsolver)
    end

    stages_efs_true = _build_stages_efs_true(S)
    stages_efs_false = _build_stages_efs_false(S)
    output = _build_output(S)
    adaptive = _build_adaptive(S)

    return quote
        $setup
        if tab.explicit_first_stage
            if integrator.f isa SplitFunction && tab.fsal && !repeat_step && !integrator.last_stepfail
                f_impl(zs[1], integrator.uprev, p, integrator.t)
                zs[1] .*= dt
            else
                @.. broadcast = false zs[1] = dt * integrator.fsalfirst
            end
            if integrator.f isa SplitFunction
                @.. broadcast = false ks[1] = dt * integrator.fsalfirst - zs[1]
            end
            $stages_efs_true
        else
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
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
            $stages_efs_false
        end

        $output

        step_limiter!(u, integrator, p, t + dt)

        $adaptive

        if integrator.f isa SplitFunction
            integrator.f(integrator.fsallast, u, p, t + dt)
        elseif tab.explicit_fsallast
            integrator.f(integrator.fsallast, u, p, t + tab.fsallast_c * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        else
            @.. broadcast = false integrator.fsallast = zs[$S] / dt
        end
    end
end

_zsym(i::Int) = Symbol("z", i)
_ksym(i::Int) = Symbol("k", i)

function _build_oop_accum_named(base, A::Symbol, zname::Function, i::Int, jrange)
    if isempty(jrange)
        return base
    end
    terms = [:($A[$i, $j] * $(zname(j))) for j in jrange]
    return Expr(:call, :+, base, terms...)
end

function _build_oop_alpha_guess_named(i::Int)
    if i == 1
        return :(zero(u))
    end
    terms = [:(α[$i][$j] * $(_zsym(j))) for j in 1:(i - 1)]
    return length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
end

function _build_oop_stages_named(S::Int; efs_true::Bool)
    body = quote end
    for i in 2:S
        zi = _zsym(i)
        ki = _ksym(i)
        z1 = _zsym(1)
        tmp_rhs = _build_oop_accum_named(:uprev, :Ai, _zsym, i, 1:(i - 1))
        ks_rhs = _build_oop_accum_named(:tmp, :Ae, _ksym, i, 1:(i - 1))
        α_rhs = _build_oop_alpha_guess_named(i)

        guess_block = if efs_true
            quote
                if integrator.f isa SplitFunction
                    z_guess = $z1
                elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
                    z_guess = tab.const_stage_guess[$i]
                elseif !isempty(α) && !iszero(α[$i])
                    z_guess = $α_rhs
                else
                    z_guess = zero(u)
                end
            end
        else
            quote
                if !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
                    z_guess = tab.const_stage_guess[$i]
                elseif !isempty(α) && !iszero(α[$i])
                    z_guess = $α_rhs
                else
                    z_guess = zero(u)
                end
            end
        end

        stage = quote
            tmp = $tmp_rhs
            if integrator.f isa SplitFunction
                tmp = $ks_rhs
            end
            $guess_block
            nlsolver.z = z_guess
            nlsolver.tmp = tmp
            nlsolver.c = c[$i]
            $zi = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            $(i < S ? quote
                if integrator.f isa SplitFunction
                    u_stage = tmp + γ * $zi
                    $ki = dt * f2(u_stage, p, t + c[$i] * dt)
                    integrator.stats.nf2 += 1
                end
            end : nothing)
        end
        append!(body.args, stage.args)
    end
    return body
end

function _build_oop_output_named(S::Int)
    u_terms = [:(bi[$i] * $(_zsym(i))) for i in 1:S]
    u_rhs = Expr(:call, :+, :uprev, u_terms...)
    u_split_terms = vcat(u_terms, [:(be[$i] * $(_ksym(i))) for i in 1:S])
    u_split_rhs = Expr(:call, :+, :uprev, u_split_terms...)
    zS = _zsym(S); kS = _ksym(S)
    return quote
        if integrator.f isa SplitFunction
            u = tmp + γ * $zS
            $kS = dt * f2(u, p, t + dt)
            integrator.stats.nf2 += 1
            u = $u_split_rhs
        elseif tab.stiffly_accurate
            u = nlsolver.tmp + γ * $zS
        else
            u = $u_rhs
        end
    end
end

function _build_oop_adaptive_named(S::Int)
    btilde_terms = [:(btilde[$i] * $(_zsym(i))) for i in 1:S]
    btilde_rhs = Expr(:call, :+, btilde_terms...)
    ebtilde_terms = [:(ebtilde[$i] * $(_ksym(i))) for i in 1:S]
    ebtilde_rhs = Expr(:call, :+, ebtilde_terms...)
    return quote
        if integrator.opts.adaptive && !isempty(btilde)
            tmp = $btilde_rhs
            if integrator.f isa SplitFunction && !isempty(ebtilde)
                tmp = tmp + $ebtilde_rhs
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
    end
end

@generated function _perform_step_oop!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{S, T, T2}
    ) where {S, T, T2}
    z1 = _zsym(1); k1 = _ksym(1); zS = _zsym(S)

    decls = quote end
    for i in 1:S
        push!(decls.args, :($(_zsym(i)) = zero(u)))
    end
    for i in 1:S
        push!(decls.args, :($(_ksym(i)) = zero(u)))
    end

    setup = quote
        (; t, dt, uprev, u, p) = integrator
        nlsolver = cache.nlsolver
        Ai = tab.Ai
        bi = tab.bi
        Ae = tab.Ae
        be = tab.be
        c = tab.c
        btilde = tab.btilde
        ebtilde = tab.ebtilde
        α = tab.α
        alg = unwrap_alg(integrator, true)
        γ = Ai[$S, $S]

        f2 = nothing
        if integrator.f isa SplitFunction
            f_impl = integrator.f.f1
            f2 = integrator.f.f2
        else
            f_impl = integrator.f
        end

        markfirststage!(nlsolver)
        $decls
    end

    stages_efs_true = _build_oop_stages_named(S; efs_true = true)
    stages_efs_false = _build_oop_stages_named(S; efs_true = false)
    output = _build_oop_output_named(S)
    adaptive = _build_oop_adaptive_named(S)

    return quote
        $setup
        if tab.explicit_first_stage
            if integrator.f isa SplitFunction
                $z1 = dt * f_impl(uprev, p, t)
            else
                $z1 = dt * integrator.fsalfirst
            end
            if integrator.f isa SplitFunction
                $k1 = dt * integrator.fsalfirst - $z1
            end
            $stages_efs_true
        else
            if integrator.success_iter > 0 && !integrator.reeval_fsal &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                    alg.extrapolant == :interpolant
                current_extrapolant!(u, t + dt, integrator)
                $z1 = u - uprev
            elseif tab.stage1_extrapolation &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm} &&
                    alg.extrapolant == :linear
                $z1 = dt * integrator.fsalfirst
            else
                $z1 = zero(u)
            end
            nlsolver.z = $z1
            nlsolver.tmp = uprev
            nlsolver.c = c[1]
            $z1 = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            $stages_efs_false
        end

        $output

        $adaptive

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
            integrator.fsallast = $zS ./ dt
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
        end
        integrator.u = u
    end
end

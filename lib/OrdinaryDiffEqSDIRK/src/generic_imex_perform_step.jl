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

# Order-limited Hermite extrapolation (power form, truncated to order q) of the previous
# step, mapped to the z-form guess z_i = (U_i - tmp)/γ. q=3 is the full Hermite.
function _build_interp_guess(i::Int, mode::Symbol)
    qsel = if mode === :max_order
        :(q_pred = 3)
    elseif mode === :cutoff_order
        :(q_pred = c[$i] <= 1 // 2 ? 3 : 1)
    else # :variable_order
        :(q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
    end
    return quote
        OrdinaryDiffEqCore.SciMLBase.addsteps!(integrator)
        if _uses_hermite_interp(alg)
            Θ_pred = (t + c[$i] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
            dtp_pred = integrator.t - integrator.tprev
            $qsel
            @.. broadcast = false zs[$i] = integrator.uprev2 +
                Θ_pred * dtp_pred * integrator.k[1]
            if q_pred >= 2
                @.. broadcast = false zs[$i] = zs[$i] +
                    Θ_pred^2 * (
                    3 * (integrator.uprev - integrator.uprev2) -
                        2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                )
            end
            if q_pred >= 3
                @.. broadcast = false zs[$i] = zs[$i] +
                    Θ_pred^3 * (
                    -2 * (integrator.uprev - integrator.uprev2) +
                        dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                )
            end
        else
            # Non-Hermite interpolant: use the full dense extrapolant (no order limiting).
            current_extrapolant!(zs[$i], t + c[$i] * dt, integrator)
        end
        @.. broadcast = false zs[$i] = (zs[$i] - tmp) * inv(γ)
    end
end

function _build_predictor_guess(i::Int, αfused)
    maxorder = _build_interp_guess(i, :max_order)
    varorder = _build_interp_guess(i, :variable_order)
    cutoff = _build_interp_guess(i, :cutoff_order)
    return quote
        if predictor == Predictor.Trivial
            fill!(zs[$i], zero(eltype(u)))
        elseif predictor == Predictor.CopyPrev && $i > 1
            copyto!(zs[$i], zs[$i - 1])
        elseif predictor == Predictor.StageExtrap && $i > 2
            @.. broadcast = false zs[$i] = zs[$i - 1] +
                (zs[$i - 1] - zs[$i - 2]) * ((c[$i] - c[$i - 1]) / (c[$i - 1] - c[$i - 2]))
        elseif predictor == Predictor.StageExtrap && $i > 1
            copyto!(zs[$i], zs[$i - 1])
        elseif predictor in
                (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) &&
                (integrator.success_iter == 0 || integrator.reeval_fsal)
            fill!(zs[$i], zero(eltype(u)))
        elseif predictor == Predictor.MaxOrder
            $maxorder
        elseif predictor == Predictor.VariableOrder
            $varorder
        elseif predictor == Predictor.CutoffOrder
            $cutoff
        elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
            fill!(zs[$i], tab.const_stage_guess[$i])
        elseif !isempty(α) && !iszero(α[$i])
            $αfused
        else
            fill!(zs[$i], zero(eltype(u)))
        end
    end
end

function _build_initial_guess_efs_false(i::Int)
    αfused = _build_alpha_guess(i)
    return _build_predictor_guess(i, αfused)
end

function _build_initial_guess_efs_true(i::Int)
    αfused = _build_alpha_guess(i)
    predbody = _build_predictor_guess(i, αfused)
    return quote
        if integrator.f isa SplitFunction && split_guess[$i] > 0
            copyto!(zs[$i], zs[split_guess[$i]])
        else
            $predbody
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
            $(
                i == 2 ? :(
                        if reuse_W_at_stage2
                            isnewton(nlsolver) && set_new_W!(nlsolver, false)
                    end
                    ) : nothing
            )
            $(
                i < S ? quote
                        if integrator.f isa SplitFunction
                            @.. broadcast = false u = tmp + γ * zs[$i]
                            f2(ks[$i], u, p, t + c[$i] * dt)
                            ks[$i] .*= dt
                            integrator.stats.nf2 += 1
                    end
                    end : nothing
            )
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

@generated function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXCache,
        tab::ESDIRKIMEXTableau{S, T, T2, :standard}, t
    ) where {S, T, T2}
    btilde_terms = [:(btilde[$i] * zs[$i]) for i in 1:S]
    btilde_rhs = Expr(:call, :+, btilde_terms...)
    ebtilde_terms = [:(ebtilde[$i] * ks[$i]) for i in 1:S]
    ebtilde_rhs = Expr(:call, :+, ebtilde_terms...)
    return quote
        (; zs, ks, atmp, nlsolver) = cache
        (; tmp) = nlsolver
        btilde = tab.btilde
        ebtilde = tab.ebtilde
        alg = unwrap_alg(integrator, true)
        uprev = integrator.uprev
        u = integrator.u
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

_build_adaptive(::Int) = :(calculate_error_estimate!(integrator, cache, tab, t))

function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXCache,
        tab::ESDIRKIMEXTableau{1, T, T2, :ie_dd2}, t
    ) where {T, T2}
    (; uprev, u, dt) = integrator
    (; atmp, zs) = cache
    return if integrator.opts.adaptive && integrator.success_iter > 0
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
        scratch = zs[1]
        @.. broadcast = false scratch = r * integrator.opts.internalnorm(
            (u - uprev) / dt1 - (uprev - uprev2) / dt2, t
        )
        calculate_residuals!(
            atmp, scratch, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    else
        OrdinaryDiffEqCore.set_EEst!(integrator, 1)
    end
end

function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXConstantCache,
        tab::ESDIRKIMEXTableau{1, T, T2, :ie_dd2}, t,
        zs::NTuple{1}, ks::NTuple{1}
    ) where {T, T2}
    (; uprev, u, dt) = integrator
    return if integrator.opts.adaptive && integrator.success_iter > 0
        uprev2 = integrator.uprev2
        tprev = integrator.tprev
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12
        r = c * dt^2
        tmp = r *
            integrator.opts.internalnorm.((u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    else
        OrdinaryDiffEqCore.set_EEst!(integrator, 1)
    end
end

@muladd function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXCache,
        tab::ESDIRKIMEXTableau{2, T, T2, :trap_dd3}, t
    ) where {T, T2}
    (; uprev, u, dt) = integrator
    (; atmp, nlsolver) = cache
    (; tmp) = nlsolver
    if integrator.opts.adaptive
        if integrator.iter > 2
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2
            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12
            r = c * dt^3 / 2
            @.. broadcast = false tmp = r * integrator.opts.internalnorm(
                (
                    (
                        (u - uprev) / dt1 - (uprev - uprev2) / dt2
                    )
                        - (
                        (uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4
                    )
                ) / dt5,
                t
            )
            calculate_residuals!(
                atmp, tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
            if OrdinaryDiffEqCore.get_EEst(integrator) <= 1
                copyto!(cache.uprev3, uprev2)
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
            copyto!(cache.uprev3, integrator.uprev2)
            cache.tprev2 = integrator.tprev
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end
end

@muladd function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXConstantCache,
        tab::ESDIRKIMEXTableau{2, T, T2, :trap_dd3}, t
    ) where {T, T2}
    (; uprev, u, dt) = integrator
    if integrator.opts.adaptive
        if integrator.iter > 2
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2
            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12
            r = c * dt^3 / 2
            DD31 = (u - uprev) / dt1 - (uprev - uprev2) / dt2
            DD30 = (uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4
            tmp = r * integrator.opts.internalnorm((DD31 - DD30) / dt5, t)
            atmp = calculate_residuals(
                tmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
            if OrdinaryDiffEqCore.get_EEst(integrator) <= 1
                cache.uprev3 = uprev2
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
            cache.uprev3 = integrator.uprev2
            cache.tprev2 = integrator.tprev
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end
end

function _build_trap_iip_body(::Int)
    return :(
        @muladd begin
            (; t, dt, uprev, u, p) = integrator
            (; atmp, nlsolver, step_limiter!) = cache
            (; z, tmp) = nlsolver
            f = integrator.f
            mass_matrix = f.mass_matrix

            γ = 1 // 2
            γdt = γ * dt
            markfirststage!(nlsolver)

            @.. broadcast = false z = uprev
            invγdt = inv(γdt)
            if mass_matrix === I
                @.. broadcast = false tmp = uprev * invγdt + integrator.fsalfirst
            else
                mul!(u, mass_matrix, uprev)
                @.. broadcast = false tmp = u * invγdt + integrator.fsalfirst
            end
            nlsolver.α = 1
            nlsolver.γ = γ
            nlsolver.method = COEFFICIENT_MULTISTEP
            z = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            @.. broadcast = false u = z

            step_limiter!(u, integrator, p, t + dt)

            calculate_error_estimate!(integrator, cache, tab, t)

            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            f(integrator.fsallast, u, p, t + dt)
        end
    )
end

function _build_trap_oop_body(::Int)
    return :(
        @muladd begin
            (; t, dt, uprev, u, p) = integrator
            nlsolver = cache.nlsolver
            f = integrator.f
            γ = 1 // 2
            γdt = γ * dt
            markfirststage!(nlsolver)

            nlsolver.z = uprev
            if f.mass_matrix === I
                nlsolver.tmp = @.. broadcast = false uprev * inv(γdt) + integrator.fsalfirst
            else
                nlsolver.tmp = (f.mass_matrix * uprev) .* inv(γdt) .+ integrator.fsalfirst
            end
            nlsolver.α = 1
            nlsolver.γ = γ
            nlsolver.method = COEFFICIENT_MULTISTEP
            u = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            integrator.u = u

            calculate_error_estimate!(integrator, cache, tab, t)

            integrator.fsallast = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            integrator.k[1] = integrator.fsalfirst
            integrator.k[2] = integrator.fsallast
        end
    )
end

@generated function _perform_step_iip!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{S, T, T2, E}
    ) where {S, T, T2, E}
    if E === :trap_dd3
        return _build_trap_iip_body(S)
    end
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
        predictor = _predictor(alg)
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
            if integrator.f isa SplitFunction && issplit(alg) && tab.fsal && !repeat_step && !integrator.last_stepfail
                f_impl(zs[1], integrator.uprev, p, integrator.t)
                zs[1] .*= dt
            else
                @.. broadcast = false zs[1] = dt * integrator.fsalfirst
            end
            if integrator.f isa SplitFunction && issplit(alg)
                @.. broadcast = false ks[1] = dt * integrator.fsalfirst - zs[1]
            end
            $stages_efs_true
        else
            if integrator.success_iter > 0 && !integrator.reeval_fsal &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm, ImplicitEuler, Trapezoid} &&
                    predictor == Predictor.MaxOrder
                current_extrapolant!(u, t + dt, integrator)
                @.. broadcast = false zs[1] = u - uprev
            elseif tab.stage1_extrapolation &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm, ImplicitEuler, Trapezoid} &&
                    predictor == Predictor.Linear
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

        if integrator.f isa SplitFunction && issplit(alg)
            integrator.f(integrator.fsallast, u, p, t + dt)
        elseif tab.explicit_fsallast
            integrator.f(integrator.fsallast, u, p, t + tab.fsallast_c * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        else
            @.. broadcast = false integrator.fsallast = zs[$S] / dt
        end

        $(
            E === :ie_dd2 ? quote
                    if integrator.opts.adaptive && integrator.differential_vars !== nothing
                        @.. broadcast = false atmp = ifelse(cache.algebraic_vars, integrator.fsallast, false) /
                        integrator.opts.abstol
                        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + integrator.opts.internalnorm(atmp, t))
                end
                end : nothing
        )
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

function _build_oop_predictor_menu(i::Int, α_rhs)
    z_prev = _zsym(i - 1)
    upred = :(
        integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
            (
            q_pred >= 2 ?
                Θ_pred^2 * (
                    3 * (integrator.uprev - integrator.uprev2) -
                    2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                ) : zero(u)
        ) +
            (
            q_pred >= 3 ?
                Θ_pred^3 * (
                    -2 * (integrator.uprev - integrator.uprev2) +
                    dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                ) : zero(u)
        )
    )
    stage_extrap = i > 2 ?
        :(
            $z_prev + ($z_prev - $(_zsym(i - 2))) *
            ((c[$i] - c[$i - 1]) / (c[$i - 1] - c[$i - 2]))
        ) : z_prev
    return quote
        if predictor == Predictor.Trivial
            z_guess = zero(u)
        elseif predictor == Predictor.CopyPrev
            z_guess = $z_prev
        elseif predictor == Predictor.StageExtrap
            z_guess = $stage_extrap
        elseif predictor in
                (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) &&
                (integrator.success_iter == 0 || integrator.reeval_fsal)
            z_guess = zero(u)
        elseif predictor in
                (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
            OrdinaryDiffEqCore.SciMLBase.addsteps!(integrator)
            Θ_pred = (t + c[$i] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
            dtp_pred = integrator.t - integrator.tprev
            q_pred = predictor == Predictor.MaxOrder ? 3 :
                predictor == Predictor.CutoffOrder ? (c[$i] <= 1 // 2 ? 3 : 1) :
                (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
            z_guess = _uses_hermite_interp(alg) ? ($upred - tmp) * inv(γ) : zero(u)
        elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[$i])
            z_guess = tab.const_stage_guess[$i]
        elseif !isempty(α) && !iszero(α[$i])
            z_guess = $α_rhs
        else
            z_guess = zero(u)
        end
    end
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
        menu = _build_oop_predictor_menu(i, α_rhs)

        guess_block = if efs_true
            quote
                if integrator.f isa SplitFunction
                    z_guess = $z1
                else
                    $menu
                end
            end
        else
            menu
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
            $(
                i < S ? quote
                        if integrator.f isa SplitFunction
                            u_stage = tmp + γ * $zi
                            $ki = dt * f2(u_stage, p, t + c[$i] * dt)
                            integrator.stats.nf2 += 1
                    end
                    end : nothing
            )
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

@generated function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXConstantCache,
        tab::ESDIRKIMEXTableau{S, T, T2, :standard}, t,
        zs::NTuple{S}, ks::NTuple{S}
    ) where {S, T, T2}
    btilde_terms = [:(btilde[$i] * zs[$i]) for i in 1:S]
    btilde_rhs = Expr(:call, :+, btilde_terms...)
    ebtilde_terms = [:(ebtilde[$i] * ks[$i]) for i in 1:S]
    ebtilde_rhs = Expr(:call, :+, ebtilde_terms...)
    return quote
        nlsolver = cache.nlsolver
        btilde = tab.btilde
        ebtilde = tab.ebtilde
        alg = unwrap_alg(integrator, true)
        uprev = integrator.uprev
        u = integrator.u
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

function _build_oop_adaptive_named(S::Int)
    zs_tuple = Expr(:tuple, (_zsym(i) for i in 1:S)...)
    ks_tuple = Expr(:tuple, (_ksym(i) for i in 1:S)...)
    return :(calculate_error_estimate!(integrator, cache, tab, t, $zs_tuple, $ks_tuple))
end

@generated function _perform_step_oop!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{S, T, T2, E}
    ) where {S, T, T2, E}
    if E === :trap_dd3
        return _build_trap_oop_body(S)
    end
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
        predictor = _predictor(alg)
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
            if integrator.f isa SplitFunction && issplit(alg)
                $z1 = dt * f_impl(uprev, p, t)
            else
                $z1 = dt * integrator.fsalfirst
            end
            if integrator.f isa SplitFunction && issplit(alg)
                $k1 = dt * integrator.fsalfirst - $z1
            end
            $stages_efs_true
        else
            if integrator.success_iter > 0 && !integrator.reeval_fsal &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm, ImplicitEuler, Trapezoid} &&
                    predictor == Predictor.MaxOrder
                current_extrapolant!(u, t + dt, integrator)
                $z1 = u - uprev
            elseif tab.stage1_extrapolation &&
                    alg isa Union{OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm, OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm, ImplicitEuler, Trapezoid} &&
                    predictor == Predictor.Linear
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

        integrator.u = u

        $adaptive

        if integrator.f isa SplitFunction && issplit(alg)
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
        $(
            E === :ie_dd2 ? quote
                    if integrator.opts.adaptive && integrator.differential_vars !== nothing
                        atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                        integrator.opts.abstol
                        OrdinaryDiffEqCore.set_EEst!(integrator, OrdinaryDiffEqCore.get_EEst(integrator) + integrator.opts.internalnorm(atmp, t))
                end
                end : nothing
        )
    end
end

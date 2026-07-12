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

# ===========================================================================
# Generic ESDIRK/IMEX perform_step bodies
# ===========================================================================
#
# The :standard and :ie_dd2 tableau flavors share the same step structure;
# they differ only in their error-estimate computation, which is inlined at
# the bottom of each body via a compile-time `E === :ie_dd2` branch (E is the
# tableau's type parameter so the branch is constant-folded).
#
# The :trap_dd3 flavor (Trapezoid with divided-difference EEst) has its own
# perform_step body further down — its stage structure differs enough that
# a hand-written specialization is clearer.
#
# Why no @generated / no Val{S}: the per-stage code below is laid out as a
# `if s >= i ... end` ladder for i ∈ 2..MAX_ESDIRKIMEX_STAGES. Each branch
# contains a fused per-stage broadcast accumulator with literal indices.
# The function compiles ONCE for every (alg type, T, T2, E) combination —
# stage count S is purely runtime — which avoids the @generated cache
# invalidation problems that motivated this rewrite.

@muladd function _perform_step_iip!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{T, T2, E}
    ) where {T, T2, E}
    (; t, dt, uprev, u, p) = integrator
    (; zs, ks, atmp, nlsolver, step_limiter!) = cache
    (; tmp) = nlsolver
    Ai = tab.Ai
    bi = tab.bi
    Ae = tab.Ae
    be = tab.be
    c = tab.c
    ce = tab.ce
    btilde = tab.btilde
    ebtilde = tab.ebtilde
    α = tab.α
    reuse_W_at_stage2 = tab.reuse_W_at_stage2
    split_guess = tab.split_guess
    alg = unwrap_alg(integrator, true)
    predictor = _predictor(alg)
    s = tab.s
    γ = Ai[s, s]

    f2 = nothing
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    # ---------------- Stage 1 ----------------
    if tab.explicit_first_stage
        if integrator.f isa SplitFunction && issplit(alg) && tab.fsal &&
                !repeat_step && !integrator.last_stepfail
            f_impl(zs[1], integrator.uprev, p, integrator.t)
            zs[1] .*= dt
        else
            @.. broadcast = false zs[1] = dt * integrator.fsalfirst
        end
        if integrator.f isa SplitFunction && issplit(alg)
            @.. broadcast = false ks[1] = dt * integrator.fsalfirst - zs[1]
        end
    else
        # Implicit first stage requires nlsolve; seed it from the requested
        # predictor. The IE tableau (E === :ie_dd2) is the only configuration
        # that reaches this branch — Trapezoid and the other Newton-SDIRKs set
        # `explicit_first_stage = true` and take the `if` arm above.
        #
        # `MaxOrder` extrapolates the previous step's interpolant, available only
        # once a step has succeeded and outside an fsal reeval; `Linear` uses
        # z = dt·f(uprev). BDF callers reuse this tableau as a first-step
        # bootstrap (e.g. ABDF2 via `cache.eulercache`) without a `predictor`
        # field of their own: `_predictor` reports `Trivial`, but
        # `!hasproperty(alg, :predictor)` keeps them on the linear bootstrap seed
        # from #3694 that preserves their convergence order under the loose
        # NonlinearSolveAlg `iter==1 && ndz<1e-5` early-exit at small dt. A
        # genuine `Trivial` request from a predictor-carrying alg falls through
        # to the zero seed.
        if E === :ie_dd2 && predictor == Predictor.MaxOrder &&
                integrator.success_iter > 0 && !integrator.reeval_fsal
            current_extrapolant!(u, t + dt, integrator)
            @.. broadcast = false zs[1] = u - uprev
        elseif E === :ie_dd2 && tab.stage1_extrapolation &&
                (predictor == Predictor.Linear || !hasproperty(alg, :predictor))
            @.. broadcast = false zs[1] = dt * integrator.fsalfirst
        else
            zs[1] .= zero(eltype(zs[1]))
        end
        nlsolver.z = zs[1]
        nlsolver.tmp = uprev
        nlsolver.c = c[1]
        zs[1] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        # All implicit stages share γ on the diagonal; reuse the W from stage 1.
        isnewton(nlsolver) && set_new_W!(nlsolver, false)
        # Non-ESDIRK IMEX (e.g. IMEX-SSP, BHR): evaluate f2 at the solved stage-1
        # value uprev + γ·zs[1]. Stage 1's explicit abscissa is t (Ae[1,:] = 0).
        if integrator.f isa SplitFunction && issplit(alg) && !isempty(ks)
            @.. broadcast = false u = uprev + γ * zs[1]
            f2(ks[1], u, p, t)
            ks[1] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    # ---------------- Stages 2..s (inlined) ----------------

    if s >= 2
        @.. broadcast = false tmp = uprev + Ai[2, 1] * zs[1]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[2, 1] * ks[1]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[2] > 0
            copyto!(zs[2], zs[split_guess[2]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[2], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[2], zs[1])
            elseif predictor == Predictor.StageExtrap
                copyto!(zs[2], zs[1])
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[2], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[2] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[2] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[2], t + c[2] * dt, integrator)
                end
                @.. broadcast = false zs[2] = (zs[2] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[2] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[2] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[2] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[2], t + c[2] * dt, integrator)
                end
                @.. broadcast = false zs[2] = (zs[2] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[2] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[2] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[2] = zs[2] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[2], t + c[2] * dt, integrator)
                end
                @.. broadcast = false zs[2] = (zs[2] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[2])
                fill!(zs[2], tab.const_stage_guess[2])
            elseif !isempty(α) && !iszero(α[2])
                @.. broadcast = false zs[2] = α[2][1] * zs[1]
            else
                fill!(zs[2], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[2]
        nlsolver.tmp = tmp
        nlsolver.c = c[2]
        zs[2] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if reuse_W_at_stage2
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
        end
        if s > 2 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[2]
            f2(ks[2], u, p, t + ce[2] * dt)
            ks[2] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 3
        @.. broadcast = false tmp = uprev + Ai[3, 1] * zs[1] + Ai[3, 2] * zs[2]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[3, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[3, 2] * ks[2]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[3] > 0
            copyto!(zs[3], zs[split_guess[3]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[3], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[3], zs[2])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[3] = zs[2] + (zs[2] - zs[1]) * ((c[3] - c[2]) / (c[2] - c[1]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[3], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[3] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[3] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[3], t + c[3] * dt, integrator)
                end
                @.. broadcast = false zs[3] = (zs[3] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[3] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[3] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[3] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[3], t + c[3] * dt, integrator)
                end
                @.. broadcast = false zs[3] = (zs[3] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[3] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[3] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[3] = zs[3] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[3], t + c[3] * dt, integrator)
                end
                @.. broadcast = false zs[3] = (zs[3] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[3])
                fill!(zs[3], tab.const_stage_guess[3])
            elseif !isempty(α) && !iszero(α[3])
                @.. broadcast = false zs[3] = α[3][1] * zs[1] + α[3][2] * zs[2]
            else
                fill!(zs[3], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[3]
        nlsolver.tmp = tmp
        nlsolver.c = c[3]
        zs[3] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 3 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[3]
            f2(ks[3], u, p, t + ce[3] * dt)
            ks[3] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 4
        @.. broadcast = false tmp = uprev + Ai[4, 1] * zs[1] + Ai[4, 2] * zs[2] + Ai[4, 3] * zs[3]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[4, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[4, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[4, 3] * ks[3]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[4] > 0
            copyto!(zs[4], zs[split_guess[4]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[4], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[4], zs[3])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[4] = zs[3] + (zs[3] - zs[2]) * ((c[4] - c[3]) / (c[3] - c[2]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[4], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[4] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[4] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[4], t + c[4] * dt, integrator)
                end
                @.. broadcast = false zs[4] = (zs[4] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[4] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[4] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[4] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[4], t + c[4] * dt, integrator)
                end
                @.. broadcast = false zs[4] = (zs[4] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[4] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[4] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[4] = zs[4] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[4], t + c[4] * dt, integrator)
                end
                @.. broadcast = false zs[4] = (zs[4] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[4])
                fill!(zs[4], tab.const_stage_guess[4])
            elseif !isempty(α) && !iszero(α[4])
                @.. broadcast = false zs[4] = α[4][1] * zs[1] + α[4][2] * zs[2] + α[4][3] * zs[3]
            else
                fill!(zs[4], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[4]
        nlsolver.tmp = tmp
        nlsolver.c = c[4]
        zs[4] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 4 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[4]
            f2(ks[4], u, p, t + ce[4] * dt)
            ks[4] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 5
        @.. broadcast = false tmp = uprev + Ai[5, 1] * zs[1] + Ai[5, 2] * zs[2] + Ai[5, 3] * zs[3] + Ai[5, 4] * zs[4]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[5, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[5, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[5, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[5, 4] * ks[4]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[5] > 0
            copyto!(zs[5], zs[split_guess[5]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[5], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[5], zs[4])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[5] = zs[4] + (zs[4] - zs[3]) * ((c[5] - c[4]) / (c[4] - c[3]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[5], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[5] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[5] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[5], t + c[5] * dt, integrator)
                end
                @.. broadcast = false zs[5] = (zs[5] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[5] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[5] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[5] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[5], t + c[5] * dt, integrator)
                end
                @.. broadcast = false zs[5] = (zs[5] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[5] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[5] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[5] = zs[5] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[5], t + c[5] * dt, integrator)
                end
                @.. broadcast = false zs[5] = (zs[5] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[5])
                fill!(zs[5], tab.const_stage_guess[5])
            elseif !isempty(α) && !iszero(α[5])
                @.. broadcast = false zs[5] = α[5][1] * zs[1] + α[5][2] * zs[2] + α[5][3] * zs[3] + α[5][4] * zs[4]
            else
                fill!(zs[5], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[5]
        nlsolver.tmp = tmp
        nlsolver.c = c[5]
        zs[5] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 5 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[5]
            f2(ks[5], u, p, t + ce[5] * dt)
            ks[5] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 6
        @.. broadcast = false tmp = uprev + Ai[6, 1] * zs[1] + Ai[6, 2] * zs[2] + Ai[6, 3] * zs[3] + Ai[6, 4] * zs[4] + Ai[6, 5] * zs[5]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[6, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[6, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[6, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[6, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[6, 5] * ks[5]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[6] > 0
            copyto!(zs[6], zs[split_guess[6]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[6], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[6], zs[5])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[6] = zs[5] + (zs[5] - zs[4]) * ((c[6] - c[5]) / (c[5] - c[4]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[6], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[6] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[6] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[6], t + c[6] * dt, integrator)
                end
                @.. broadcast = false zs[6] = (zs[6] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[6] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[6] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[6] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[6], t + c[6] * dt, integrator)
                end
                @.. broadcast = false zs[6] = (zs[6] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[6] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[6] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[6] = zs[6] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[6], t + c[6] * dt, integrator)
                end
                @.. broadcast = false zs[6] = (zs[6] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[6])
                fill!(zs[6], tab.const_stage_guess[6])
            elseif !isempty(α) && !iszero(α[6])
                @.. broadcast = false zs[6] = α[6][1] * zs[1] + α[6][2] * zs[2] + α[6][3] * zs[3] + α[6][4] * zs[4] + α[6][5] * zs[5]
            else
                fill!(zs[6], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[6]
        nlsolver.tmp = tmp
        nlsolver.c = c[6]
        zs[6] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 6 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[6]
            f2(ks[6], u, p, t + ce[6] * dt)
            ks[6] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 7
        @.. broadcast = false tmp = uprev + Ai[7, 1] * zs[1] + Ai[7, 2] * zs[2] + Ai[7, 3] * zs[3] + Ai[7, 4] * zs[4] + Ai[7, 5] * zs[5] + Ai[7, 6] * zs[6]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[7, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[7, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[7, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[7, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[7, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[7, 6] * ks[6]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[7] > 0
            copyto!(zs[7], zs[split_guess[7]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[7], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[7], zs[6])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[7] = zs[6] + (zs[6] - zs[5]) * ((c[7] - c[6]) / (c[6] - c[5]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[7], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[7] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[7] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[7], t + c[7] * dt, integrator)
                end
                @.. broadcast = false zs[7] = (zs[7] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[7] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[7] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[7] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[7], t + c[7] * dt, integrator)
                end
                @.. broadcast = false zs[7] = (zs[7] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[7] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[7] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[7] = zs[7] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[7], t + c[7] * dt, integrator)
                end
                @.. broadcast = false zs[7] = (zs[7] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[7])
                fill!(zs[7], tab.const_stage_guess[7])
            elseif !isempty(α) && !iszero(α[7])
                @.. broadcast = false zs[7] = α[7][1] * zs[1] + α[7][2] * zs[2] + α[7][3] * zs[3] + α[7][4] * zs[4] + α[7][5] * zs[5] + α[7][6] * zs[6]
            else
                fill!(zs[7], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[7]
        nlsolver.tmp = tmp
        nlsolver.c = c[7]
        zs[7] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 7 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[7]
            f2(ks[7], u, p, t + ce[7] * dt)
            ks[7] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 8
        @.. broadcast = false tmp = uprev + Ai[8, 1] * zs[1] + Ai[8, 2] * zs[2] + Ai[8, 3] * zs[3] + Ai[8, 4] * zs[4] + Ai[8, 5] * zs[5] + Ai[8, 6] * zs[6] + Ai[8, 7] * zs[7]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[8, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[8, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[8, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[8, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[8, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[8, 6] * ks[6]
            @.. broadcast = false tmp = tmp + Ae[8, 7] * ks[7]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[8] > 0
            copyto!(zs[8], zs[split_guess[8]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[8], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[8], zs[7])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[8] = zs[7] + (zs[7] - zs[6]) * ((c[8] - c[7]) / (c[7] - c[6]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[8], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[8] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[8] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[8], t + c[8] * dt, integrator)
                end
                @.. broadcast = false zs[8] = (zs[8] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[8] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[8] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[8] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[8], t + c[8] * dt, integrator)
                end
                @.. broadcast = false zs[8] = (zs[8] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[8] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[8] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[8] = zs[8] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[8], t + c[8] * dt, integrator)
                end
                @.. broadcast = false zs[8] = (zs[8] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[8])
                fill!(zs[8], tab.const_stage_guess[8])
            elseif !isempty(α) && !iszero(α[8])
                @.. broadcast = false zs[8] = α[8][1] * zs[1] + α[8][2] * zs[2] + α[8][3] * zs[3] + α[8][4] * zs[4] + α[8][5] * zs[5] + α[8][6] * zs[6] + α[8][7] * zs[7]
            else
                fill!(zs[8], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[8]
        nlsolver.tmp = tmp
        nlsolver.c = c[8]
        zs[8] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 8 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[8]
            f2(ks[8], u, p, t + ce[8] * dt)
            ks[8] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 9
        @.. broadcast = false tmp = uprev + Ai[9, 1] * zs[1] + Ai[9, 2] * zs[2] + Ai[9, 3] * zs[3] + Ai[9, 4] * zs[4] + Ai[9, 5] * zs[5] + Ai[9, 6] * zs[6] + Ai[9, 7] * zs[7] + Ai[9, 8] * zs[8]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[9, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[9, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[9, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[9, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[9, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[9, 6] * ks[6]
            @.. broadcast = false tmp = tmp + Ae[9, 7] * ks[7]
            @.. broadcast = false tmp = tmp + Ae[9, 8] * ks[8]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[9] > 0
            copyto!(zs[9], zs[split_guess[9]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[9], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[9], zs[8])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[9] = zs[8] + (zs[8] - zs[7]) * ((c[9] - c[8]) / (c[8] - c[7]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[9], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[9] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[9] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[9], t + c[9] * dt, integrator)
                end
                @.. broadcast = false zs[9] = (zs[9] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[9] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[9] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[9] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[9], t + c[9] * dt, integrator)
                end
                @.. broadcast = false zs[9] = (zs[9] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[9] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[9] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[9] = zs[9] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[9], t + c[9] * dt, integrator)
                end
                @.. broadcast = false zs[9] = (zs[9] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[9])
                fill!(zs[9], tab.const_stage_guess[9])
            elseif !isempty(α) && !iszero(α[9])
                @.. broadcast = false zs[9] = α[9][1] * zs[1] + α[9][2] * zs[2] + α[9][3] * zs[3] + α[9][4] * zs[4] + α[9][5] * zs[5] + α[9][6] * zs[6] + α[9][7] * zs[7] + α[9][8] * zs[8]
            else
                fill!(zs[9], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[9]
        nlsolver.tmp = tmp
        nlsolver.c = c[9]
        zs[9] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 9 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[9]
            f2(ks[9], u, p, t + ce[9] * dt)
            ks[9] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 10
        @.. broadcast = false tmp = uprev + Ai[10, 1] * zs[1] + Ai[10, 2] * zs[2] + Ai[10, 3] * zs[3] + Ai[10, 4] * zs[4] + Ai[10, 5] * zs[5] + Ai[10, 6] * zs[6] + Ai[10, 7] * zs[7] + Ai[10, 8] * zs[8] + Ai[10, 9] * zs[9]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[10, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[10, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[10, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[10, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[10, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[10, 6] * ks[6]
            @.. broadcast = false tmp = tmp + Ae[10, 7] * ks[7]
            @.. broadcast = false tmp = tmp + Ae[10, 8] * ks[8]
            @.. broadcast = false tmp = tmp + Ae[10, 9] * ks[9]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[10] > 0
            copyto!(zs[10], zs[split_guess[10]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[10], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[10], zs[9])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[10] = zs[9] + (zs[9] - zs[8]) * ((c[10] - c[9]) / (c[9] - c[8]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[10], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[10] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[10] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[10], t + c[10] * dt, integrator)
                end
                @.. broadcast = false zs[10] = (zs[10] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[10] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[10] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[10] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[10], t + c[10] * dt, integrator)
                end
                @.. broadcast = false zs[10] = (zs[10] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[10] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[10] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[10] = zs[10] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[10], t + c[10] * dt, integrator)
                end
                @.. broadcast = false zs[10] = (zs[10] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[10])
                fill!(zs[10], tab.const_stage_guess[10])
            elseif !isempty(α) && !iszero(α[10])
                @.. broadcast = false zs[10] = α[10][1] * zs[1] + α[10][2] * zs[2] + α[10][3] * zs[3] + α[10][4] * zs[4] + α[10][5] * zs[5] + α[10][6] * zs[6] + α[10][7] * zs[7] + α[10][8] * zs[8] + α[10][9] * zs[9]
            else
                fill!(zs[10], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[10]
        nlsolver.tmp = tmp
        nlsolver.c = c[10]
        zs[10] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 10 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[10]
            f2(ks[10], u, p, t + ce[10] * dt)
            ks[10] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 11
        @.. broadcast = false tmp = uprev + Ai[11, 1] * zs[1] + Ai[11, 2] * zs[2] + Ai[11, 3] * zs[3] + Ai[11, 4] * zs[4] + Ai[11, 5] * zs[5] + Ai[11, 6] * zs[6] + Ai[11, 7] * zs[7] + Ai[11, 8] * zs[8] + Ai[11, 9] * zs[9] + Ai[11, 10] * zs[10]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[11, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[11, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[11, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[11, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[11, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[11, 6] * ks[6]
            @.. broadcast = false tmp = tmp + Ae[11, 7] * ks[7]
            @.. broadcast = false tmp = tmp + Ae[11, 8] * ks[8]
            @.. broadcast = false tmp = tmp + Ae[11, 9] * ks[9]
            @.. broadcast = false tmp = tmp + Ae[11, 10] * ks[10]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[11] > 0
            copyto!(zs[11], zs[split_guess[11]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[11], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[11], zs[10])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[11] = zs[10] + (zs[10] - zs[9]) * ((c[11] - c[10]) / (c[10] - c[9]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[11], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[11] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[11] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[11], t + c[11] * dt, integrator)
                end
                @.. broadcast = false zs[11] = (zs[11] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[11] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[11] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[11] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[11], t + c[11] * dt, integrator)
                end
                @.. broadcast = false zs[11] = (zs[11] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[11] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[11] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[11] = zs[11] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[11], t + c[11] * dt, integrator)
                end
                @.. broadcast = false zs[11] = (zs[11] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[11])
                fill!(zs[11], tab.const_stage_guess[11])
            elseif !isempty(α) && !iszero(α[11])
                @.. broadcast = false zs[11] = α[11][1] * zs[1] + α[11][2] * zs[2] + α[11][3] * zs[3] + α[11][4] * zs[4] + α[11][5] * zs[5] + α[11][6] * zs[6] + α[11][7] * zs[7] + α[11][8] * zs[8] + α[11][9] * zs[9] + α[11][10] * zs[10]
            else
                fill!(zs[11], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[11]
        nlsolver.tmp = tmp
        nlsolver.c = c[11]
        zs[11] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 11 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[11]
            f2(ks[11], u, p, t + ce[11] * dt)
            ks[11] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    if s >= 12
        @.. broadcast = false tmp = uprev + Ai[12, 1] * zs[1] + Ai[12, 2] * zs[2] + Ai[12, 3] * zs[3] + Ai[12, 4] * zs[4] + Ai[12, 5] * zs[5] + Ai[12, 6] * zs[6] + Ai[12, 7] * zs[7] + Ai[12, 8] * zs[8] + Ai[12, 9] * zs[9] + Ai[12, 10] * zs[10] + Ai[12, 11] * zs[11]
        if integrator.f isa SplitFunction
            @.. broadcast = false tmp = tmp + Ae[12, 1] * ks[1]
            @.. broadcast = false tmp = tmp + Ae[12, 2] * ks[2]
            @.. broadcast = false tmp = tmp + Ae[12, 3] * ks[3]
            @.. broadcast = false tmp = tmp + Ae[12, 4] * ks[4]
            @.. broadcast = false tmp = tmp + Ae[12, 5] * ks[5]
            @.. broadcast = false tmp = tmp + Ae[12, 6] * ks[6]
            @.. broadcast = false tmp = tmp + Ae[12, 7] * ks[7]
            @.. broadcast = false tmp = tmp + Ae[12, 8] * ks[8]
            @.. broadcast = false tmp = tmp + Ae[12, 9] * ks[9]
            @.. broadcast = false tmp = tmp + Ae[12, 10] * ks[10]
            @.. broadcast = false tmp = tmp + Ae[12, 11] * ks[11]
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction && split_guess[12] > 0
            copyto!(zs[12], zs[split_guess[12]])
        else
            if predictor == Predictor.Trivial
                fill!(zs[12], zero(eltype(u)))
            elseif predictor == Predictor.CopyPrev
                copyto!(zs[12], zs[11])
            elseif predictor == Predictor.StageExtrap
                @.. broadcast = false zs[12] = zs[11] + (zs[11] - zs[10]) * ((c[12] - c[11]) / (c[11] - c[10]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                fill!(zs[12], zero(eltype(u)))
            elseif predictor == Predictor.MaxOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[12] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = 3
                    @.. broadcast = false zs[12] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[12], t + c[12] * dt, integrator)
                end
                @.. broadcast = false zs[12] = (zs[12] - tmp) * inv(γ)
            elseif predictor == Predictor.CutoffOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[12] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = c[12] <= 1 // 2 ? 3 : 1
                    @.. broadcast = false zs[12] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[12], t + c[12] * dt, integrator)
                end
                @.. broadcast = false zs[12] = (zs[12] - tmp) * inv(γ)
            elseif predictor == Predictor.VariableOrder
                SciMLBase.addsteps!(integrator)
                if _uses_hermite_interp(alg)
                    Θ_pred = (t + c[12] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                    dtp_pred = integrator.t - integrator.tprev
                    q_pred = Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1)
                    @.. broadcast = false zs[12] = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1]
                    if q_pred >= 2
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^2 * (3 * (integrator.uprev - integrator.uprev2) - 2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2])
                    end
                    if q_pred >= 3
                        @.. broadcast = false zs[12] = zs[12] + Θ_pred^3 * (-2 * (integrator.uprev - integrator.uprev2) + dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2])
                    end
                else
                    current_extrapolant!(zs[12], t + c[12] * dt, integrator)
                end
                @.. broadcast = false zs[12] = (zs[12] - tmp) * inv(γ)
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[12])
                fill!(zs[12], tab.const_stage_guess[12])
            elseif !isempty(α) && !iszero(α[12])
                @.. broadcast = false zs[12] = α[12][1] * zs[1] + α[12][2] * zs[2] + α[12][3] * zs[3] + α[12][4] * zs[4] + α[12][5] * zs[5] + α[12][6] * zs[6] + α[12][7] * zs[7] + α[12][8] * zs[8] + α[12][9] * zs[9] + α[12][10] * zs[10] + α[12][11] * zs[11]
            else
                fill!(zs[12], zero(eltype(u)))
            end
        end
        nlsolver.z = zs[12]
        nlsolver.tmp = tmp
        nlsolver.c = c[12]
        zs[12] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 12 && integrator.f isa SplitFunction
            @.. broadcast = false u = tmp + γ * zs[12]
            f2(ks[12], u, p, t + ce[12] * dt)
            ks[12] .*= dt
            integrator.stats.nf2 += 1
        end
    end


    # ---------------- Output u ----------------
    if integrator.f isa SplitFunction
        @.. broadcast = false u = tmp + γ * zs[s]
        f2(ks[s], u, p, t + dt)
        ks[s] .*= dt
        integrator.stats.nf2 += 1
        if s == 1
            @.. broadcast = false u = uprev + bi[1] * zs[1] + be[1] * ks[1]
        elseif s == 2
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + be[1] * ks[1] + be[2] * ks[2]
        elseif s == 3
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3]
        elseif s == 4
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4]
        elseif s == 5
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5]
        elseif s == 6
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6]
        elseif s == 7
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7]
        elseif s == 8
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7] + be[8] * ks[8]
        elseif s == 9
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7] + be[8] * ks[8] + be[9] * ks[9]
        elseif s == 10
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7] + be[8] * ks[8] + be[9] * ks[9] + be[10] * ks[10]
        elseif s == 11
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10] + bi[11] * zs[11] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7] + be[8] * ks[8] + be[9] * ks[9] + be[10] * ks[10] + be[11] * ks[11]
        elseif s == 12
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10] + bi[11] * zs[11] + bi[12] * zs[12] + be[1] * ks[1] + be[2] * ks[2] + be[3] * ks[3] + be[4] * ks[4] + be[5] * ks[5] + be[6] * ks[6] + be[7] * ks[7] + be[8] * ks[8] + be[9] * ks[9] + be[10] * ks[10] + be[11] * ks[11] + be[12] * ks[12]
        end
    elseif tab.stiffly_accurate
        @.. broadcast = false u = tmp + γ * zs[s]
    else
        if s == 1
            @.. broadcast = false u = uprev + bi[1] * zs[1]
        elseif s == 2
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2]
        elseif s == 3
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3]
        elseif s == 4
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4]
        elseif s == 5
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5]
        elseif s == 6
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6]
        elseif s == 7
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7]
        elseif s == 8
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8]
        elseif s == 9
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9]
        elseif s == 10
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10]
        elseif s == 11
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10] + bi[11] * zs[11]
        elseif s == 12
            @.. broadcast = false u = uprev + bi[1] * zs[1] + bi[2] * zs[2] + bi[3] * zs[3] + bi[4] * zs[4] + bi[5] * zs[5] + bi[6] * zs[6] + bi[7] * zs[7] + bi[8] * zs[8] + bi[9] * zs[9] + bi[10] * zs[10] + bi[11] * zs[11] + bi[12] * zs[12]
        end
    end

    step_limiter!(u, integrator, p, t + dt)

    # ---------------- Error estimate ----------------
    if E === :standard
        if integrator.opts.adaptive && !isempty(btilde)
            if s == 1
                @.. broadcast = false tmp = btilde[1] * zs[1]
            elseif s == 2
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2]
            elseif s == 3
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3]
            elseif s == 4
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4]
            elseif s == 5
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5]
            elseif s == 6
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6]
            elseif s == 7
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7]
            elseif s == 8
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7] + btilde[8] * zs[8]
            elseif s == 9
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7] + btilde[8] * zs[8] + btilde[9] * zs[9]
            elseif s == 10
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7] + btilde[8] * zs[8] + btilde[9] * zs[9] + btilde[10] * zs[10]
            elseif s == 11
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7] + btilde[8] * zs[8] + btilde[9] * zs[9] + btilde[10] * zs[10] + btilde[11] * zs[11]
            elseif s == 12
                @.. broadcast = false tmp = btilde[1] * zs[1] + btilde[2] * zs[2] + btilde[3] * zs[3] + btilde[4] * zs[4] + btilde[5] * zs[5] + btilde[6] * zs[6] + btilde[7] * zs[7] + btilde[8] * zs[8] + btilde[9] * zs[9] + btilde[10] * zs[10] + btilde[11] * zs[11] + btilde[12] * zs[12]
            end
            if integrator.f isa SplitFunction && !isempty(ebtilde)
                for j in 1:s
                    @.. broadcast = false tmp = tmp + ebtilde[j] * ks[j]
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
    elseif E === :ie_dd2
        # Divided-difference-of-order-2 estimate for ImplicitEuler (s == 1).
        if integrator.opts.adaptive && integrator.success_iter > 0
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            cc = 7 / 12
            r = cc * dt^2
            scratch = zs[1]
            @.. broadcast = false scratch = r * integrator.opts.internalnorm(
                (u - uprev) / dt1 - (uprev - uprev2) / dt2, t
            )
            calculate_residuals!(
                atmp, scratch, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        elseif integrator.opts.adaptive
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end

    # ---------------- fsallast ----------------
    if integrator.f isa SplitFunction && issplit(alg)
        integrator.f(integrator.fsallast, u, p, t + dt)
    elseif tab.explicit_fsallast
        integrator.f(integrator.fsallast, u, p, t + tab.fsallast_c * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    else
        @.. broadcast = false integrator.fsallast = zs[s] / dt
    end

    # ---------------- :ie_dd2-specific DAE EEst tail ----------------
    if E === :ie_dd2
        if integrator.opts.adaptive && integrator.differential_vars !== nothing
            @.. broadcast = false atmp =
                ifelse(cache.algebraic_vars, integrator.fsallast, false) /
                integrator.opts.abstol
            OrdinaryDiffEqCore.set_EEst!(
                integrator,
                OrdinaryDiffEqCore.get_EEst(integrator) +
                    integrator.opts.internalnorm(atmp, t)
            )
        end
    end
    return nothing
end


@muladd function _perform_step_oop!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{T, T2, E}
    ) where {T, T2, E}
    (; t, dt, uprev, u, p) = integrator
    nlsolver = cache.nlsolver
    Ai = tab.Ai
    bi = tab.bi
    Ae = tab.Ae
    be = tab.be
    c = tab.c
    ce = tab.ce
    btilde = tab.btilde
    ebtilde = tab.ebtilde
    α = tab.α
    reuse_W_at_stage2 = tab.reuse_W_at_stage2
    alg = unwrap_alg(integrator, true)
    predictor = _predictor(alg)
    s = tab.s
    γ = Ai[s, s]

    f2 = nothing
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    # Named locals z1..z<MAX>, k1..k<MAX>. Each z_i is filled by its matching
    # `if s >= i` block, and only referenced inside the matching `if s == k`
    # output ladder, so the `local` declaration is enough — no allocation
    # cost for slots beyond runtime `s`.
    #
    # The k_i story is messier: when `integrator.f isa SplitFunction` is true
    # but `issplit(alg)` is false (e.g. Kvaerno4 fed a SplitODEProblem), the
    # stage-1 setup never assigns k1, yet stages 2..s still enter the
    # `tmp += Ae[i,j]*kj` accumulator (Ae is all zeros in that case, but Julia
    # still needs k1 defined to evaluate `Ae[2,1] * k1`). Pre-allocating all
    # k_i to `zero(u)` whenever f is a SplitFunction matches master's
    # @generated behavior and costs at most 12 array zeros per step for the
    # split path; for the common non-split path the `local` declaration
    # leaves them un-allocated.
    local z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12
    local k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12
    if integrator.f isa SplitFunction
        k1 = zero(u); k2 = zero(u); k3 = zero(u); k4 = zero(u)
        k5 = zero(u); k6 = zero(u); k7 = zero(u); k8 = zero(u)
        k9 = zero(u); k10 = zero(u); k11 = zero(u); k12 = zero(u)
    end
    tmp = uprev

    # ---------------- Stage 1 ----------------
    if tab.explicit_first_stage
        if integrator.f isa SplitFunction && issplit(alg)
            z1 = dt * f_impl(uprev, p, t)
            k1 = dt * integrator.fsalfirst - z1
        else
            z1 = dt * integrator.fsalfirst
        end
    else
        # See the matching branch in `_perform_step_iip!` above. `u` is immutable
        # here (scalar / SVector out-of-place), so the MaxOrder seed uses the
        # allocating `current_extrapolant` rather than the in-place variant.
        if E === :ie_dd2 && predictor == Predictor.MaxOrder &&
                integrator.success_iter > 0 && !integrator.reeval_fsal
            z1 = current_extrapolant(t + dt, integrator) - uprev
        elseif E === :ie_dd2 && tab.stage1_extrapolation &&
                (predictor == Predictor.Linear || !hasproperty(alg, :predictor))
            z1 = dt * integrator.fsalfirst
        else
            z1 = zero(u)
        end
        nlsolver.z = z1
        nlsolver.tmp = uprev
        nlsolver.c = c[1]
        z1 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        # Non-ESDIRK IMEX (e.g. IMEX-SSP, BHR): evaluate f2 at the solved stage-1
        # value uprev + γ·z1. Stage 1's explicit abscissa is t (Ae[1,:] = 0).
        if integrator.f isa SplitFunction && issplit(alg)
            u_stage1 = uprev + γ * z1
            k1 = dt * f2(u_stage1, p, t)
            integrator.stats.nf2 += 1
        end
    end

    # ---------------- Stages 2..s (inlined) ----------------

    if s >= 2
        tmp = uprev + Ai[2, 1] * z1
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[2, 1] * k1
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z1
            elseif predictor == Predictor.StageExtrap
                z_guess = z1
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[2] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[2] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[2])
                z_guess = tab.const_stage_guess[2]
            elseif !isempty(α) && !iszero(α[2])
                z_guess = α[2][1] * z1
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[2]
        z2 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if reuse_W_at_stage2
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
        end
        if s > 2 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z2
            k2 = dt * f2(u_stage, p, t + ce[2] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 3
        tmp = uprev + Ai[3, 1] * z1 + Ai[3, 2] * z2
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[3, 1] * k1 + Ae[3, 2] * k2
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z2
            elseif predictor == Predictor.StageExtrap
                z_guess = z2 + (z2 - z1) * ((c[3] - c[2]) / (c[2] - c[1]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[3] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[3] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[3])
                z_guess = tab.const_stage_guess[3]
            elseif !isempty(α) && !iszero(α[3])
                z_guess = α[3][1] * z1 + α[3][2] * z2
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[3]
        z3 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 3 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z3
            k3 = dt * f2(u_stage, p, t + ce[3] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 4
        tmp = uprev + Ai[4, 1] * z1 + Ai[4, 2] * z2 + Ai[4, 3] * z3
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[4, 1] * k1 + Ae[4, 2] * k2 + Ae[4, 3] * k3
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z3
            elseif predictor == Predictor.StageExtrap
                z_guess = z3 + (z3 - z2) * ((c[4] - c[3]) / (c[3] - c[2]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[4] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[4] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[4])
                z_guess = tab.const_stage_guess[4]
            elseif !isempty(α) && !iszero(α[4])
                z_guess = α[4][1] * z1 + α[4][2] * z2 + α[4][3] * z3
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[4]
        z4 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 4 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z4
            k4 = dt * f2(u_stage, p, t + ce[4] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 5
        tmp = uprev + Ai[5, 1] * z1 + Ai[5, 2] * z2 + Ai[5, 3] * z3 + Ai[5, 4] * z4
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[5, 1] * k1 + Ae[5, 2] * k2 + Ae[5, 3] * k3 + Ae[5, 4] * k4
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z4
            elseif predictor == Predictor.StageExtrap
                z_guess = z4 + (z4 - z3) * ((c[5] - c[4]) / (c[4] - c[3]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[5] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[5] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[5])
                z_guess = tab.const_stage_guess[5]
            elseif !isempty(α) && !iszero(α[5])
                z_guess = α[5][1] * z1 + α[5][2] * z2 + α[5][3] * z3 + α[5][4] * z4
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[5]
        z5 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 5 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z5
            k5 = dt * f2(u_stage, p, t + ce[5] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 6
        tmp = uprev + Ai[6, 1] * z1 + Ai[6, 2] * z2 + Ai[6, 3] * z3 + Ai[6, 4] * z4 + Ai[6, 5] * z5
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[6, 1] * k1 + Ae[6, 2] * k2 + Ae[6, 3] * k3 + Ae[6, 4] * k4 + Ae[6, 5] * k5
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z5
            elseif predictor == Predictor.StageExtrap
                z_guess = z5 + (z5 - z4) * ((c[6] - c[5]) / (c[5] - c[4]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[6] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[6] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[6])
                z_guess = tab.const_stage_guess[6]
            elseif !isempty(α) && !iszero(α[6])
                z_guess = α[6][1] * z1 + α[6][2] * z2 + α[6][3] * z3 + α[6][4] * z4 + α[6][5] * z5
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[6]
        z6 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 6 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z6
            k6 = dt * f2(u_stage, p, t + ce[6] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 7
        tmp = uprev + Ai[7, 1] * z1 + Ai[7, 2] * z2 + Ai[7, 3] * z3 + Ai[7, 4] * z4 + Ai[7, 5] * z5 + Ai[7, 6] * z6
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[7, 1] * k1 + Ae[7, 2] * k2 + Ae[7, 3] * k3 + Ae[7, 4] * k4 + Ae[7, 5] * k5 + Ae[7, 6] * k6
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z6
            elseif predictor == Predictor.StageExtrap
                z_guess = z6 + (z6 - z5) * ((c[7] - c[6]) / (c[6] - c[5]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[7] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[7] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[7])
                z_guess = tab.const_stage_guess[7]
            elseif !isempty(α) && !iszero(α[7])
                z_guess = α[7][1] * z1 + α[7][2] * z2 + α[7][3] * z3 + α[7][4] * z4 + α[7][5] * z5 + α[7][6] * z6
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[7]
        z7 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 7 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z7
            k7 = dt * f2(u_stage, p, t + ce[7] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 8
        tmp = uprev + Ai[8, 1] * z1 + Ai[8, 2] * z2 + Ai[8, 3] * z3 + Ai[8, 4] * z4 + Ai[8, 5] * z5 + Ai[8, 6] * z6 + Ai[8, 7] * z7
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[8, 1] * k1 + Ae[8, 2] * k2 + Ae[8, 3] * k3 + Ae[8, 4] * k4 + Ae[8, 5] * k5 + Ae[8, 6] * k6 + Ae[8, 7] * k7
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z7
            elseif predictor == Predictor.StageExtrap
                z_guess = z7 + (z7 - z6) * ((c[8] - c[7]) / (c[7] - c[6]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[8] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[8] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[8])
                z_guess = tab.const_stage_guess[8]
            elseif !isempty(α) && !iszero(α[8])
                z_guess = α[8][1] * z1 + α[8][2] * z2 + α[8][3] * z3 + α[8][4] * z4 + α[8][5] * z5 + α[8][6] * z6 + α[8][7] * z7
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[8]
        z8 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 8 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z8
            k8 = dt * f2(u_stage, p, t + ce[8] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 9
        tmp = uprev + Ai[9, 1] * z1 + Ai[9, 2] * z2 + Ai[9, 3] * z3 + Ai[9, 4] * z4 + Ai[9, 5] * z5 + Ai[9, 6] * z6 + Ai[9, 7] * z7 + Ai[9, 8] * z8
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[9, 1] * k1 + Ae[9, 2] * k2 + Ae[9, 3] * k3 + Ae[9, 4] * k4 + Ae[9, 5] * k5 + Ae[9, 6] * k6 + Ae[9, 7] * k7 + Ae[9, 8] * k8
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z8
            elseif predictor == Predictor.StageExtrap
                z_guess = z8 + (z8 - z7) * ((c[9] - c[8]) / (c[8] - c[7]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[9] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[9] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[9])
                z_guess = tab.const_stage_guess[9]
            elseif !isempty(α) && !iszero(α[9])
                z_guess = α[9][1] * z1 + α[9][2] * z2 + α[9][3] * z3 + α[9][4] * z4 + α[9][5] * z5 + α[9][6] * z6 + α[9][7] * z7 + α[9][8] * z8
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[9]
        z9 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 9 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z9
            k9 = dt * f2(u_stage, p, t + ce[9] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 10
        tmp = uprev + Ai[10, 1] * z1 + Ai[10, 2] * z2 + Ai[10, 3] * z3 + Ai[10, 4] * z4 + Ai[10, 5] * z5 + Ai[10, 6] * z6 + Ai[10, 7] * z7 + Ai[10, 8] * z8 + Ai[10, 9] * z9
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[10, 1] * k1 + Ae[10, 2] * k2 + Ae[10, 3] * k3 + Ae[10, 4] * k4 + Ae[10, 5] * k5 + Ae[10, 6] * k6 + Ae[10, 7] * k7 + Ae[10, 8] * k8 + Ae[10, 9] * k9
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z9
            elseif predictor == Predictor.StageExtrap
                z_guess = z9 + (z9 - z8) * ((c[10] - c[9]) / (c[9] - c[8]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[10] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[10] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[10])
                z_guess = tab.const_stage_guess[10]
            elseif !isempty(α) && !iszero(α[10])
                z_guess = α[10][1] * z1 + α[10][2] * z2 + α[10][3] * z3 + α[10][4] * z4 + α[10][5] * z5 + α[10][6] * z6 + α[10][7] * z7 + α[10][8] * z8 + α[10][9] * z9
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[10]
        z10 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 10 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z10
            k10 = dt * f2(u_stage, p, t + ce[10] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 11
        tmp = uprev + Ai[11, 1] * z1 + Ai[11, 2] * z2 + Ai[11, 3] * z3 + Ai[11, 4] * z4 + Ai[11, 5] * z5 + Ai[11, 6] * z6 + Ai[11, 7] * z7 + Ai[11, 8] * z8 + Ai[11, 9] * z9 + Ai[11, 10] * z10
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[11, 1] * k1 + Ae[11, 2] * k2 + Ae[11, 3] * k3 + Ae[11, 4] * k4 + Ae[11, 5] * k5 + Ae[11, 6] * k6 + Ae[11, 7] * k7 + Ae[11, 8] * k8 + Ae[11, 9] * k9 + Ae[11, 10] * k10
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z10
            elseif predictor == Predictor.StageExtrap
                z_guess = z10 + (z10 - z9) * ((c[11] - c[10]) / (c[10] - c[9]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[11] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[11] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[11])
                z_guess = tab.const_stage_guess[11]
            elseif !isempty(α) && !iszero(α[11])
                z_guess = α[11][1] * z1 + α[11][2] * z2 + α[11][3] * z3 + α[11][4] * z4 + α[11][5] * z5 + α[11][6] * z6 + α[11][7] * z7 + α[11][8] * z8 + α[11][9] * z9 + α[11][10] * z10
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[11]
        z11 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 11 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z11
            k11 = dt * f2(u_stage, p, t + ce[11] * dt)
            integrator.stats.nf2 += 1
        end
    end

    if s >= 12
        tmp = uprev + Ai[12, 1] * z1 + Ai[12, 2] * z2 + Ai[12, 3] * z3 + Ai[12, 4] * z4 + Ai[12, 5] * z5 + Ai[12, 6] * z6 + Ai[12, 7] * z7 + Ai[12, 8] * z8 + Ai[12, 9] * z9 + Ai[12, 10] * z10 + Ai[12, 11] * z11
        if integrator.f isa SplitFunction
            tmp = tmp + Ae[12, 1] * k1 + Ae[12, 2] * k2 + Ae[12, 3] * k3 + Ae[12, 4] * k4 + Ae[12, 5] * k5 + Ae[12, 6] * k6 + Ae[12, 7] * k7 + Ae[12, 8] * k8 + Ae[12, 9] * k9 + Ae[12, 10] * k10 + Ae[12, 11] * k11
        end
        if tab.explicit_first_stage && integrator.f isa SplitFunction
            z_guess = z1
        else
            if predictor == Predictor.Trivial
                z_guess = zero(u)
            elseif predictor == Predictor.CopyPrev
                z_guess = z11
            elseif predictor == Predictor.StageExtrap
                z_guess = z11 + (z11 - z10) * ((c[12] - c[11]) / (c[11] - c[10]))
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder) && (integrator.success_iter == 0 || integrator.reeval_fsal)
                z_guess = zero(u)
            elseif predictor in (Predictor.MaxOrder, Predictor.VariableOrder, Predictor.CutoffOrder)
                SciMLBase.addsteps!(integrator)
                Θ_pred = (t + c[12] * dt - integrator.tprev) / (integrator.t - integrator.tprev)
                dtp_pred = integrator.t - integrator.tprev
                q_pred = predictor == Predictor.MaxOrder ? 3 :
                    predictor == Predictor.CutoffOrder ? (c[12] <= 1 // 2 ? 3 : 1) :
                    (Θ_pred <= 3 // 2 ? 3 : (Θ_pred <= 5 // 2 ? 2 : 1))
                if _uses_hermite_interp(alg)
                    u_pred = integrator.uprev2 + Θ_pred * dtp_pred * integrator.k[1] +
                        (
                        q_pred >= 2 ?
                            Θ_pred^2 * (
                                3 * (integrator.uprev - integrator.uprev2) -
                                2 * dtp_pred * integrator.k[1] - dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    ) +
                        (
                        q_pred >= 3 ?
                            Θ_pred^3 * (
                                -2 * (integrator.uprev - integrator.uprev2) +
                                dtp_pred * integrator.k[1] + dtp_pred * integrator.k[2]
                            ) :
                            zero(u)
                    )
                    z_guess = (u_pred - tmp) * inv(γ)
                else
                    z_guess = zero(u)
                end
            elseif !isempty(tab.const_stage_guess) && !iszero(tab.const_stage_guess[12])
                z_guess = tab.const_stage_guess[12]
            elseif !isempty(α) && !iszero(α[12])
                z_guess = α[12][1] * z1 + α[12][2] * z2 + α[12][3] * z3 + α[12][4] * z4 + α[12][5] * z5 + α[12][6] * z6 + α[12][7] * z7 + α[12][8] * z8 + α[12][9] * z9 + α[12][10] * z10 + α[12][11] * z11
            else
                z_guess = zero(u)
            end
        end
        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[12]
        z12 = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if s > 12 && integrator.f isa SplitFunction
            u_stage = tmp + γ * z12
            k12 = dt * f2(u_stage, p, t + ce[12] * dt)
            integrator.stats.nf2 += 1
        end
    end


    # ---------------- Output u ----------------
    if integrator.f isa SplitFunction
        if s == 1
            u_last = tmp + γ * z1
            k1 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + be[1] * k1
        elseif s == 2
            u_last = tmp + γ * z2
            k2 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + be[1] * k1 + be[2] * k2
        elseif s == 3
            u_last = tmp + γ * z3
            k3 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + be[1] * k1 + be[2] * k2 + be[3] * k3
        elseif s == 4
            u_last = tmp + γ * z4
            k4 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4
        elseif s == 5
            u_last = tmp + γ * z5
            k5 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5
        elseif s == 6
            u_last = tmp + γ * z6
            k6 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6
        elseif s == 7
            u_last = tmp + γ * z7
            k7 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7
        elseif s == 8
            u_last = tmp + γ * z8
            k8 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7 + be[8] * k8
        elseif s == 9
            u_last = tmp + γ * z9
            k9 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7 + be[8] * k8 + be[9] * k9
        elseif s == 10
            u_last = tmp + γ * z10
            k10 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7 + be[8] * k8 + be[9] * k9 + be[10] * k10
        elseif s == 11
            u_last = tmp + γ * z11
            k11 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10 + bi[11] * z11 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7 + be[8] * k8 + be[9] * k9 + be[10] * k10 + be[11] * k11
        elseif s == 12
            u_last = tmp + γ * z12
            k12 = dt * f2(u_last, p, t + dt)
            integrator.stats.nf2 += 1
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10 + bi[11] * z11 + bi[12] * z12 + be[1] * k1 + be[2] * k2 + be[3] * k3 + be[4] * k4 + be[5] * k5 + be[6] * k6 + be[7] * k7 + be[8] * k8 + be[9] * k9 + be[10] * k10 + be[11] * k11 + be[12] * k12
        end
    elseif tab.stiffly_accurate
        if s == 1
            u = nlsolver.tmp + γ * z1
        elseif s == 2
            u = nlsolver.tmp + γ * z2
        elseif s == 3
            u = nlsolver.tmp + γ * z3
        elseif s == 4
            u = nlsolver.tmp + γ * z4
        elseif s == 5
            u = nlsolver.tmp + γ * z5
        elseif s == 6
            u = nlsolver.tmp + γ * z6
        elseif s == 7
            u = nlsolver.tmp + γ * z7
        elseif s == 8
            u = nlsolver.tmp + γ * z8
        elseif s == 9
            u = nlsolver.tmp + γ * z9
        elseif s == 10
            u = nlsolver.tmp + γ * z10
        elseif s == 11
            u = nlsolver.tmp + γ * z11
        elseif s == 12
            u = nlsolver.tmp + γ * z12
        end
    else
        if s == 1
            u = uprev + bi[1] * z1
        elseif s == 2
            u = uprev + bi[1] * z1 + bi[2] * z2
        elseif s == 3
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3
        elseif s == 4
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4
        elseif s == 5
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5
        elseif s == 6
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6
        elseif s == 7
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7
        elseif s == 8
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8
        elseif s == 9
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9
        elseif s == 10
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10
        elseif s == 11
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10 + bi[11] * z11
        elseif s == 12
            u = uprev + bi[1] * z1 + bi[2] * z2 + bi[3] * z3 + bi[4] * z4 + bi[5] * z5 + bi[6] * z6 + bi[7] * z7 + bi[8] * z8 + bi[9] * z9 + bi[10] * z10 + bi[11] * z11 + bi[12] * z12
        end
    end

    integrator.u = u

    # ---------------- Error estimate ----------------
    if E === :standard
        if integrator.opts.adaptive && !isempty(btilde)
            local tmp_est
            if s == 1
                tmp_est = btilde[1] * z1
            elseif s == 2
                tmp_est = btilde[1] * z1 + btilde[2] * z2
            elseif s == 3
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3
            elseif s == 4
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4
            elseif s == 5
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5
            elseif s == 6
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6
            elseif s == 7
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7
            elseif s == 8
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7 + btilde[8] * z8
            elseif s == 9
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7 + btilde[8] * z8 + btilde[9] * z9
            elseif s == 10
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7 + btilde[8] * z8 + btilde[9] * z9 + btilde[10] * z10
            elseif s == 11
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7 + btilde[8] * z8 + btilde[9] * z9 + btilde[10] * z10 + btilde[11] * z11
            elseif s == 12
                tmp_est = btilde[1] * z1 + btilde[2] * z2 + btilde[3] * z3 + btilde[4] * z4 + btilde[5] * z5 + btilde[6] * z6 + btilde[7] * z7 + btilde[8] * z8 + btilde[9] * z9 + btilde[10] * z10 + btilde[11] * z11 + btilde[12] * z12
            end
            if integrator.f isa SplitFunction && !isempty(ebtilde)
                if s == 1
                    tmp_est = tmp_est + ebtilde[1] * k1
                elseif s == 2
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2
                elseif s == 3
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3
                elseif s == 4
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4
                elseif s == 5
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5
                elseif s == 6
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6
                elseif s == 7
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7
                elseif s == 8
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7 + ebtilde[8] * k8
                elseif s == 9
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7 + ebtilde[8] * k8 + ebtilde[9] * k9
                elseif s == 10
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7 + ebtilde[8] * k8 + ebtilde[9] * k9 + ebtilde[10] * k10
                elseif s == 11
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7 + ebtilde[8] * k8 + ebtilde[9] * k9 + ebtilde[10] * k10 + ebtilde[11] * k11
                elseif s == 12
                    tmp_est = tmp_est + ebtilde[1] * k1 + ebtilde[2] * k2 + ebtilde[3] * k3 + ebtilde[4] * k4 + ebtilde[5] * k5 + ebtilde[6] * k6 + ebtilde[7] * k7 + ebtilde[8] * k8 + ebtilde[9] * k9 + ebtilde[10] * k10 + ebtilde[11] * k11 + ebtilde[12] * k12
                end
            end
            if isnewton(nlsolver) && _esdirk_smooth_est(alg)
                integrator.stats.nsolve += 1
                est = _reshape(get_W(nlsolver) \ _vec(tmp_est), axes(tmp_est))
            else
                est = tmp_est
            end
            atmp = calculate_residuals(
                est, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        end
    elseif E === :ie_dd2
        # Divided-difference-of-order-2 estimate for ImplicitEuler (s == 1).
        if integrator.opts.adaptive && integrator.success_iter > 0
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            cc = 7 / 12
            r = cc * dt^2
            tmp_est = r * integrator.opts.internalnorm.(
                (u - uprev) / dt1 - (uprev - uprev2) / dt2, t
            )
            atmp = calculate_residuals(
                tmp_est, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        elseif integrator.opts.adaptive
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    end

    # ---------------- fsallast + k1/k2 ----------------
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
        if s == 1
            integrator.fsallast = z1 ./ dt
        elseif s == 2
            integrator.fsallast = z2 ./ dt
        elseif s == 3
            integrator.fsallast = z3 ./ dt
        elseif s == 4
            integrator.fsallast = z4 ./ dt
        elseif s == 5
            integrator.fsallast = z5 ./ dt
        elseif s == 6
            integrator.fsallast = z6 ./ dt
        elseif s == 7
            integrator.fsallast = z7 ./ dt
        elseif s == 8
            integrator.fsallast = z8 ./ dt
        elseif s == 9
            integrator.fsallast = z9 ./ dt
        elseif s == 10
            integrator.fsallast = z10 ./ dt
        elseif s == 11
            integrator.fsallast = z11 ./ dt
        elseif s == 12
            integrator.fsallast = z12 ./ dt
        end
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end

    # ---------------- :ie_dd2-specific DAE EEst tail ----------------
    if E === :ie_dd2
        if integrator.opts.adaptive && integrator.differential_vars !== nothing
            atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) /
                integrator.opts.abstol
            OrdinaryDiffEqCore.set_EEst!(
                integrator,
                OrdinaryDiffEqCore.get_EEst(integrator) +
                    integrator.opts.internalnorm(atmp, t)
            )
        end
    end
    return nothing
end


# ===========================================================================
# :ie_dd2 — error-estimate helpers used only as references in the body above
# (the EEst lives inside _perform_step_iip! / _perform_step_oop!). Nothing
# to add here; the divided-difference-2 fallback for `success_iter == 0`
# is handled by the standard :ie_dd2 EEst block emitted by the bodies above
# returning the default OrdinaryDiffEqCore set_EEst path when relevant.
# ===========================================================================

# ===========================================================================
# :trap_dd3 — Trapezoid with divided-difference order-3 error estimate.
# Trapezoid has fixed γ = 1//2 and a stage structure compact enough that a
# hand-written specialization is clearer than going through the generic
# perform_step bodies above.
# ===========================================================================
@muladd function calculate_error_estimate!(
        integrator, cache::ESDIRKIMEXCache,
        tab::ESDIRKIMEXTableau{T, T2, :trap_dd3}, t
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
        tab::ESDIRKIMEXTableau{T, T2, :trap_dd3}, t
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

@muladd function _perform_step_iip!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{T, T2, :trap_dd3}
    ) where {T, T2}
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
    return f(integrator.fsallast, u, p, t + dt)
end

@muladd function _perform_step_oop!(
        integrator, cache, repeat_step, tab::ESDIRKIMEXTableau{T, T2, :trap_dd3}
    ) where {T, T2}
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
    return integrator.k[2]
end

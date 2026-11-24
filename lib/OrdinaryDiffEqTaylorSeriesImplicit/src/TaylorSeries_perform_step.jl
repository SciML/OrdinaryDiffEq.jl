function initialize!(integrator, cache::ImplicitTaylorConstantCache)
    reset_newton_state!(cache.stage_cache, integrator.t)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitTaylorCache)
    reset_newton_state!(cache, integrator.t)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# Also called by reinit!: clear all history, factorizations and feedback.
function reset_newton_state!(cache, t)
    cache.t = cache.tprev2 = cache.t_hist = t
    cache.jac_valid = cache.lu_valid = cache.filter_valid = false
    cache.embedded_filtered = false
    cache.newton_failed = cache.newton_limit = cache.recovery_since_accept = false
    cache.newton_rate_measured = cache.newton_slow = false
    cache.dt_cap = zero(cache.dt_cap)
    cache.cap_probe, cache.cap_streak, cache.cap_raised = 1, 0, false
    cache.ηold = cache.ηleak = one(cache.ηold)
    cache.status = VerySlowConvergence
    cache.iter = 0
    cache.ncoefjac = cache.ncheckjac = cache.nendpoint_refresh = 0
    cache.nrestart = cache.nvalidation_fail = cache.nfilter_recovery = cache.nleak = 0
    cache.reduced_check_h = zero(cache.reduced_check_h)
    fill!(cache.uprev3, zero(eltype(cache.uprev3)))
    fill!(cache.uhist, zero(eltype(cache.uhist)))
    cache.factored === nothing || reset_factored!(cache.factored)
    cache.stage_jac_valid = cache.lu_stages_valid = cache.stages_valid = false
    cache.jprev_valid = false
    cache.use_stages = cache.stage_system && !cache.stage_auto
    cache.stage_trial = cache.stage_recovery = cache.stage_handover = false
    cache.stage_probe, cache.stage_streak = 1, 0
    return nothing
end

function perform_step!(integrator, cache::ImplicitTaylorConstantCache, repeat_step = false)
    work = cache.stage_cache
    if integrator.uprev isa Number
        work.uprev[1] = integrator.uprev
        work.uprev2[1] = integrator.uprev2
        work.fsalfirst[1] = integrator.fsalfirst
    else
        copyto!(work.uprev, vec(integrator.uprev))
        copyto!(work.uprev2, vec(integrator.uprev2))
        copyto!(work.fsalfirst, vec(integrator.fsalfirst))
    end
    implicit_taylor_step!(integrator, work, integrator.t, integrator.dt, work.uprev, work.u) || return
    u = taylor_result(integrator.uprev, work.u)
    integrator.u = u
    integrator.fsallast = integrator.f(u, integrator.p, integrator.t + integrator.dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function perform_step!(integrator, cache::ImplicitTaylorCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    implicit_taylor_step!(integrator, cache, t, dt, uprev, u) || return
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

# One nonlinear solve, then the error estimate; `false` if Newton failed.
function implicit_taylor_step!(integrator, cache, t, dt, uprev, u)
    try
        implicit_taylor_stage!(integrator, cache, t, dt, uprev, u) || return false
        if integrator.opts.adaptive
            implicit_taylor_error!(integrator, cache, t, dt, uprev, u)
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
        end
    catch err
        err isa TaylorLinearSolveFailure || rethrow()
        # A failed filter can occur after the candidate has been written to u.
        # Discard that candidate and all factors that the attempt could have changed.
        copyto!(u, uprev)
        cache.jac_valid = cache.lu_valid = cache.filter_valid = false
        cache.stage_jac_valid = cache.lu_stages_valid = cache.stages_valid = false
        cache.jprev_valid = cache.embedded_filtered = false
        cache.reduced_check_h = zero(cache.reduced_check_h)
        cache.factored === nothing ||
            (cache.factored.distance[] = convert(eltype(cache.factored.correction), NaN))
        cache.status = Divergence
        cache.newton_failed = cache.recovery_since_accept = true
        integrator.opts.adaptive && record_newton_failure!(cache, dt)
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        OrdinaryDiffEqCore.set_EEst!(integrator, Inf)
        return false
    end
    return true
end

# Predict z = u_stage - uprev. History extrapolation is quadratic with three
# states, and linear near startup or when the step grows by more than 3×.
function taylor_initial_guess!(integrator, cache, alg, t, dt, uprev, pade)
    (; utilde, rhs, predictor, fsalfirst, uprev2, uprev3) = cache
    P, Q = _unwrap_val(alg.order), _unwrap_val(alg.order_q)
    order = P + Q
    # Distance from uprev to the implicit point.
    distance = pade ? dt : cache.μ * dt
    if alg.extrapolant === :extrapolate && integrator.success_iter > 0 &&
            t != integrator.tprev
        h1 = t - integrator.tprev
        h2 = integrator.tprev - cache.tprev2
        if integrator.success_iter > 1 && h2 != zero(h2) && abs(distance) <= 3 * abs(h1)
            # Quadratic through (tprev2, uprev3), (tprev, uprev2), (t, uprev)
            # in Newton form, evaluated at t + distance.
            curvature = distance * (distance + h1) / (h1 + h2)
            c1 = distance / h1 + curvature / h1
            c2 = curvature / h2
            @. utilde = c1 * (uprev - uprev2) - c2 * (uprev2 - uprev3)
        else
            scale = distance / h1
            @. utilde = scale * (uprev - uprev2)
        end
    elseif alg.extrapolant === :polynomial && pade && P > 0 && Q > 0 && order > 2
        # α A_inc + β B_inc at (uprev, t) equals h y' + h²/2 y'' for every (P, Q).
        α = (order * (order + Q - 2)) / (P * (order - 2))
        β = order / (order - 2)
        @. utilde = α * rhs + β * predictor
    elseif alg.extrapolant === :linear
        @. utilde = distance * fsalfirst
    else
        fill!(utilde, zero(eltype(utilde)))
    end
    return nothing
end

# View the solver's vector buffer in the shape of the state.
@inline shaped_like(v, u) = size(v) == size(u) ? v : _reshape(v, size(u))

# M = I + Σ_k w_k d^k ∂c_k/∂u.
function assemble_newton_matrix!(J, Jk, weights, distance)
    fill!(J, zero(eltype(J)))
    n = size(J, 1)
    @inbounds for k in 1:size(Jk, 3)
        w = weights[k + 1] * distance^k
        iszero(w) && continue
        for j in 1:n, i in 1:n
            J[i, j] += w * Jk[i, j, k]
        end
    end
    @inbounds for i in 1:n
        J[i, i] += one(eltype(J))
    end
    return J
end

# Solve z + B_inc(uprev + z) = A_inc(uprev) for Taylor-Padé, or
# z + T_inc(uprev + z, t + μh, -μh) = 0 for μ-Taylor.
# Adaptive solves require η‖Δz‖ < κ and a filtered residual check. A slow
# component can escape the correction norm when contraction rates differ.
function implicit_taylor_stage!(integrator, cache, t, dt, uprev, u)
    (; μ, κ, tmp, tmp_embedded, atmp, linsolve, J, Jk, uintermediate, rhs, rhs_embedded, predictor, utilde) = cache
    alg = unwrap_alg(integrator, true)
    (; internalnorm, abstol, reltol, adaptive) = integrator.opts
    (; maxiters) = alg
    P, Q = _unwrap_val(alg.order), _unwrap_val(alg.order_q)
    pade = is_taylor_pade(alg)
    stats = integrator.stats
    cache.embedded_filtered = false

    # Expansion center, distance and weights of the implicit polynomial.
    if pade
        tc = t + dt
        distance = dt
        weights = pade_increment_coefficients(Val(P), Val(Q))[2]
    else
        tc = t + μ * dt
        distance = -μ * dt
        weights = taylor_increment_coefficients(Val(P))
    end

    # Explicit side at (uprev, t); the copy gives the kernels their element type.
    uintermediate .= uprev
    if pade && P > 0
        cache.polynomial_prev(rhs, rhs_embedded, predictor, uintermediate, t, dt)
    else
        fill!(rhs, zero(eltype(rhs)))
        fill!(rhs_embedded, zero(eltype(rhs_embedded)))
    end

    # New step (not a retry): shift the history.
    if cache.t != t
        copyto!(cache.uprev3, cache.uhist)
        copyto!(cache.uhist, cache.uprev2)
        cache.tprev2 = cache.t_hist
        cache.t_hist = integrator.tprev
    end
    # Constant Jacobians survive accepted steps. Variable Jacobians persist
    # only through error retries. A discontinuity invalidates both.
    if integrator.derivative_discontinuity
        cache.jac_valid = cache.lu_valid = cache.filter_valid = false
        cache.stage_jac_valid = cache.lu_stages_valid = cache.stages_valid = false
        cache.jprev_valid = false
        cache.reduced_check_h = zero(cache.reduced_check_h)
        cache.ηleak = one(cache.ηleak)
        cache.factored === nothing || (cache.factored.distance[] = convert(eltype(cache.factored.correction), NaN))
    end
    if cache.stage_system
        if cache.use_stages
            new_step = cache.t != t
            cache.t = t
            return derivative_stage!(integrator, cache, t, dt, uprev, u, tc, distance, weights, new_step)
        end
        # stages = :auto on the increment form: restored if this solve hands its step over.
        cache.stage_recovery = cache.recovery_since_accept
    end
    fresh = !(adaptive && cache.jac_valid && (cache.t == t || cache.constant_jacobian))
    cache.t = t
    # Fixed steps use the assembled matrix.
    if adaptive && cache.factored !== nothing && cache.factored.active[]
        return factored_stage!(integrator, cache, cache.factored, t, dt, uprev, u, tc, distance, weights)
    end
    # Fixed steps refactorize at the predictor. Without a second-order predictor,
    # orders above two refresh the matrix again at the second correction.
    force_factorize = !adaptive
    polynomial_start = alg.extrapolant === :polynomial && pade && P > 0 && Q > 0 && P + Q > 2
    refresh_second = !adaptive && !polynomial_start && alg_order(alg) > 2
    θ_slow = convert(typeof(κ), 1 // 2)
    κ_second = κ * convert(typeof(κ), 1 // 10)
    tiny = κ * convert(typeof(κ), 1 // 1000)
    # Solution buffer, never rebound to the (not fully inferred) solution object.
    dw = linsolve.u
    correction = shaped_like(dw, utilde)
    η = ndw = ndwprev = θ = one(κ)
    converged = false
    jac_fresh = false  # the coefficient Jacobians are valid at this attempt's predictor
    predictor_jacobian_valid = cache.constant_jacobian && cache.jac_valid
    factorized = false # M was factorized in this attempt
    measured = false   # θ was measured in the last attempt
    iter = 0
    total_iter = 0
    endpoint_refresh_used = false
    validation_failed = false
    cache_retry = false
    # Restart at the predictor with a fresh linearization if reuse contracts poorly.
    for attempt in 1:2
        taylor_initial_guess!(integrator, cache, alg, t, dt, uprev, pade)
        @. uintermediate = uprev + utilde
        cache.status = VerySlowConvergence
        η = one(η)
        jac_iter = 1 # iteration at which M was assembled in this attempt
        θ_from = 3
        iter = 0
        pending_matrix = false
        jac_fresh = factorized = false
        while iter < maxiters
            iter += 1
            total_iter += 1
            stats.nnonliniter += 1
            factorize = false
            residual_ready = pending_matrix
            if pending_matrix
                # Endpoint refresh: residual and matrix at the current iterate.
                copyto!(cache.Mkept, cache.Mtmp)
                matrix_row_norms!(cache.Mrow_norm, cache.Mkept)
                copyto!(J, cache.Mtmp)
                cache.lu_valid = true
                factorize = factorized = jac_fresh = true
                pending_matrix = false
                cache.nendpoint_refresh += 1
                cache.recovery_since_accept = true
                jac_iter = iter
                θ_from = iter + 2
            elseif iter == 1 || (iter == 2 && refresh_second)
                if fresh || iter == 2
                    # Values and coefficient Jacobians share one jet evaluation.
                    cache.polynomial_jac(tmp, tmp_embedded, Jk, uintermediate, tc, distance)
                    stats.njacs += 1
                    cache.ncoefjac += 1
                    cache.jac_valid = true
                    jac_fresh = true
                    predictor_jacobian_valid = iter == 1
                else
                    cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, distance)
                    # A retry at the same predictor needs a fresh LU, not a
                    # second evaluation of its already-current derivatives.
                    jac_fresh = predictor_jacobian_valid || cache.constant_jacobian
                end
                assemble_newton_matrix!(cache.Mtmp, Jk, weights, distance)
                if force_factorize || newton_matrix_outdated(cache, matrix_cutoff(alg))
                    copyto!(cache.Mkept, cache.Mtmp)
                    matrix_row_norms!(cache.Mrow_norm, cache.Mkept)
                    copyto!(J, cache.Mtmp)
                    cache.lu_valid = true
                    factorize = factorized = true
                end
                jac_iter = iter
                θ_from = iter + 2
            else
                cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, distance)
            end
            if !residual_ready
                @. tmp = utilde + tmp - rhs
            end
            factorize && (stats.nw += 1)
            if cache.direct_lu
                if factorize && !lu_factor!(J, cache.pivJ)
                    cache.lu_valid = false
                    cache.status = Divergence
                    break
                end
                copyto!(dw, _vec(tmp))
                lu_solve!(J, cache.pivJ, dw)
                if !all(isfinite, dw)
                    cache.lu_valid = false
                    cache.status = Divergence
                    break
                end
            elseif factorize
                linres = checked_dolinsolve(integrator, linsolve; A = J, b = _vec(tmp), linu = dw)
                cache.detached_solution && copyto!(dw, linres.u)
            else
                linres = checked_dolinsolve(integrator, linsolve; b = _vec(tmp), linu = dw)
                cache.detached_solution && copyto!(dw, linres.u)
            end
            stats.nsolve += 1

            # Norm of the correction, scaled by the integration tolerances.
            iter > 1 && (ndwprev = ndw)
            @.. broadcast = false atmp = calculate_residuals(
                abs(correction), uprev, uintermediate, abstol, reltol, internalnorm, t
            )
            ndw = internalnorm(atmp, t)
            if !isfinite(ndw)
                cache.status = Divergence
                break
            end
            # If the retained LU contracts poorly, restart before applying
            # this correction. The provisional rate is discarded.
            if adaptive && jac_fresh && !factorized && iter == jac_iter + 1 &&
                    ndw >= κ_second && ndw > θ_slow * ndwprev &&
                    cache.Mtmp != cache.Mkept
                cache.status = VerySlowConvergence
                break
            end
            @. utilde = utilde - correction
            @. uintermediate = uprev + utilde

            bad_convergence = false
            if iter >= θ_from
                θ = ndw / ndwprev
                bad_convergence = θ > 2 || (
                    adaptive &&
                        (θ >= one(θ) || (!(jac_fresh && factorized) && θ > θ_slow))
                )
                η = θ / (1 - θ)
            end

            candidate = ndw <= tiny || (iter == jac_iter + 1 && ndw < κ_second) ||
                (iter >= θ_from && η >= zero(η) && η * ndw < κ)
            # Tiny corrections may stagnate at roundoff; they still go through the check.
            if bad_convergence && ndw > tiny
                cache.status = Divergence
                break
            end
            if adaptive && candidate
                degree = size(Jk, 3)
                # A nonconstant degree-one endpoint check overwrites Jk;
                # higher degrees use the separate Jcheck view.
                degree == 1 && (predictor_jacobian_valid = false)
                defect, changed, paired = check_stage_defect!(
                    integrator, cache, t, tc, distance, uprev, weights
                )
                if !isfinite(defect)
                    cache.nvalidation_fail += 1
                    validation_failed = true
                    cache.status = Divergence
                    break
                end
                if !changed && defect <= κ
                    converged = true
                    cache.embedded_filtered = paired
                    break
                end
                cache.nvalidation_fail += 1
                if endpoint_refresh_used || iter == maxiters
                    # Fallback after a failed endpoint refresh; degree-one filters
                    # share the Newton matrix and give no independent direction.
                    if endpoint_refresh_used && degree > 1 && iter < maxiters
                        converged, trials = recover_filtered_stage!(
                            integrator, cache, t, tc, distance, uprev, weights,
                            defect, maxiters - iter
                        )
                        iter += trials
                        total_iter += trials
                        # Not a chord correction: reset η.
                        η = one(η)
                        θ_from = iter + 2
                    end
                    validation_failed = !converged
                    converged || (cache.status = VerySlowConvergence)
                    break
                end
                # At most one endpoint refresh per stage solve.
                endpoint_refresh_used = true
                if degree > 1
                    predictor_jacobian_valid = false
                    cache.polynomial_jac(tmp, tmp_embedded, Jk, uintermediate, tc, distance)
                    stats.njacs += 1
                    cache.ncoefjac += 1
                    cache.jac_valid = true
                    assemble_newton_matrix!(cache.Mtmp, Jk, weights, distance)
                    @. tmp = utilde + tmp - rhs
                end
                pending_matrix = true
                η = one(η)
                continue
            elseif candidate
                converged = true
                break
            end
            if bad_convergence
                cache.status = Divergence
                break
            end
            # The remaining budget cannot reach the tolerance at this θ.
            if adaptive && iter >= θ_from && θ < one(θ) &&
                    ndw * θ^(maxiters - iter) > κ * (1 - θ)
                cache.status = VerySlowConvergence
                break
            end
        end
        measured = iter >= θ_from
        # Retry once with fresh Jacobians and factorization unless both were fresh.
        (converged || validation_failed || (jac_fresh && factorized) || attempt == 2) && break
        cache.nrestart += 1
        # The fresh attempt distinguishes factor reuse from step-size difficulty.
        cache_retry = true
        fresh = !predictor_jacobian_valid
        force_factorize = true
    end
    cache.iter = total_iter
    if !converged
        # The retry must not reuse this linearization.
        cache.jac_valid = false
        if cache.stage_system
            # :auto retries the same step in the derivative-stage form.
            stats.nnonlinconvfail += 1
            hand_over_to_stages!(cache, adaptive)
            return derivative_stage!(integrator, cache, t, dt, uprev, u, tc, distance, weights, true)
        end
        cache.newton_failed = true
        adaptive && record_newton_failure!(cache, dt)
        check_rounding_limit!(integrator, cache, alg, adaptive, dw)
        integrator.force_stepfail = true
        stats.nnonlinconvfail += 1
        return false
    end
    cache.status = Convergence
    # Only fast convergence permits growth. A validated fast retry with fresh
    # factors clears a reuse failure; other recovery persists until acceptance.
    fast = !measured || zero(η) <= η <= fast_convergence_cutoff(alg)
    cache_retry && !fast && (cache.recovery_since_accept = true)
    record_newton_convergence!(cache, measured, η, alg)
    cache.newton_limit && check_rounding_limit!(integrator, cache, alg, adaptive, dw)
    # The derivative stages' last converged stages no longer precede the next step.
    cache.stages_valid = false

    # Map the stage value to the end of the step.
    if pade
        u .= uintermediate
    else
        # μ-Taylor: forward increment from the center; the same evaluation gives
        # c_P d^P for the error estimate (at d = -h when μ = 1).
        if μ != one(μ) || (adaptive && uses_last_term(alg))
            d_end = μ == one(μ) ? -oftype(tc, dt) : (1 - μ) * dt
            cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, d_end)
        end
        if μ == one(μ)
            u .= uintermediate
        elseif alg.real_function
            @. u = real(uintermediate + tmp)
        else
            @. u = uintermediate + tmp
        end
    end
    return true
end

# Filtered residual of a converged candidate, independent of the Newton LU.
# For degree one the filter is M, so the matrix is checked for change as well.
@inline function check_stage_defect!(integrator, cache, t, tc, distance, uprev, weights)
    (; tmp, tmp_embedded, Jk, uintermediate, utilde, rhs, rhs_embedded, atmp) = cache
    (; abstol, reltol, internalnorm) = integrator.opts
    degree = size(Jk, 3)
    if cache.constant_jacobian && cache.jac_valid
        # The symbolic proof applies to the derivatives, not the polynomial
        # values: time-dependent forcing must still be evaluated at this point.
        cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, distance)
        if degree == 1
            assemble_newton_matrix!(cache.Mtmp, Jk, weights, distance)
        else
            copyto!(cache.Mtmp, @view(Jk[:, :, 1]))
        end
    else
        if degree == 1
            cache.polynomial_jac(tmp, tmp_embedded, Jk, uintermediate, tc, distance)
            assemble_newton_matrix!(cache.Mtmp, Jk, weights, distance)
        else
            # Endpoint f' into Jcheck (aliases Mtmp).
            cache.polynomial_check(tmp, tmp_embedded, cache.Jcheck, uintermediate, tc, distance)
        end
        integrator.stats.njacs += 1
        cache.ncheckjac += 1
    end
    alg = unwrap_alg(integrator, true)
    changed = degree == 1 && newton_matrix_outdated(cache, matrix_cutoff(alg))
    @. tmp = utilde + tmp - rhs
    copyto!(cache.stage_error, tmp)
    paired = degree > 1 && cache.linsolve_filter.b isa AbstractMatrix
    if paired
        @. tmp_embedded = utilde + tmp_embedded - rhs_embedded
    end
    filter_error!(
        integrator, cache, cache.stage_error, degree, weights[2], distance,
        degree > 1 ? cache.Mtmp : nothing, paired ? tmp_embedded : nothing
    )
    @.. broadcast = false atmp = calculate_residuals(
        abs(cache.stage_error), uprev, uintermediate, abstol, reltol, internalnorm, t
    )
    return internalnorm(atmp, t), changed, paired
end

# Last resort after validation failed even with a refreshed matrix: backtrack
# along the filtered residual with fixed filter factors; only a fresh
# validation accepts the result.
@noinline function recover_filtered_stage!(
        integrator, cache, t, tc, distance, uprev, weights, defect, budget
    )
    budget > 0 || return false, 0
    (; utilde, uintermediate, tmp, tmp_embedded, stage_error, atmp, κ) = cache
    (; abstol, reltol, internalnorm) = integrator.opts
    degree = size(cache.Jk, 3)
    # Rollback storage, allocated only on this path.
    start_increment = copy(utilde)
    start_state = copy(uintermediate)
    direction = copy(stage_error)
    cache.nfilter_recovery += 1
    cache.recovery_since_accept = true
    cache.embedded_filtered = false
    α = one(κ)
    trials = 0
    while trials < budget
        trials += 1
        integrator.stats.nnonliniter += 1
        @. utilde = start_increment - α * direction
        @. uintermediate = uprev + utilde
        if all(isfinite, uintermediate)
            cache.polynomial_iter(tmp, tmp_embedded, uintermediate, tc, distance)
            @. stage_error = utilde + tmp - cache.rhs
            apply_filter_factors!(integrator, cache, stage_error, degree)
            @.. broadcast = false atmp = calculate_residuals(
                abs(stage_error), uprev, start_state, abstol, reltol, internalnorm, t
            )
            trial_defect = internalnorm(atmp, t)
            # Sufficient decrease; stop at the first such trial.
            if isfinite(trial_defect) && trial_defect <= (1 - α / 2) * defect
                if trial_defect <= κ
                    checked, changed, paired = check_stage_defect!(
                        integrator, cache, t, tc, distance, uprev, weights
                    )
                    if isfinite(checked) && !changed && checked <= κ
                        cache.embedded_filtered = paired
                        # An error retry must not reuse the old linearization.
                        cache.jac_valid = cache.lu_valid = false
                        return true, trials
                    end
                    cache.nvalidation_fail += 1
                end
                break
            end
        end
        α /= 2
    end
    copyto!(utilde, start_increment)
    copyto!(uintermediate, start_state)
    # Buffers describe a rejected trial point: invalidate them.
    cache.jac_valid = cache.lu_valid = cache.filter_valid = false
    cache.embedded_filtered = false
    return false, trials
end

# Compare each row against cutoff × max(1, sum(abs, kept_row)).
# Comparing rows separately keeps the stiff rows from masking changes elsewhere.
function newton_matrix_outdated(cache, cutoff)
    (; Mkept, Mtmp) = cache
    cache.lu_valid || return true
    return rows_changed(Mtmp, Mkept, cutoff, cache.Mrow_norm, cache.row_change)
end

# Column-major scans in blocks of 32 rows, with early exit per block.
function matrix_row_norms!(norms, A)
    isempty(norms) && return nothing
    fill!(norms, zero(eltype(norms)))
    @inbounds for j in axes(A, 2)
        @simd for i in axes(A, 1)
            norms[i] += abs(A[i, j])
        end
    end
    return nothing
end

function rows_changed(new, kept, cutoff, norms, changes, scale = nothing, leak = nothing)
    isempty(norms) && return rows_changed_scan(new, kept, cutoff, scale, leak)
    n = size(kept, 1)
    fill!(changes, zero(eltype(changes)))
    @inbounds for firstrow in 1:32:n
        lastrow = min(firstrow + 31, n)
        for j in 1:n
            @simd for i in firstrow:lastrow
                value = scale === nothing ? new[i, j] : scale * new[i, j]
                changes[i] += abs(value - kept[i, j])
            end
        end
        for i in firstrow:lastrow
            changes[i] > row_change_limit(cutoff, norms[i], leak) && return true
        end
    end
    return false
end

# Without stored row norms: the kept row's norm is formed in the same pass.
function rows_changed_scan(new, kept, cutoff, scale, leak)
    n = size(new, 1)
    @inbounds for i in 1:n
        row_change = row_kept = zero(real(eltype(new)))
        for j in 1:n
            k = kept[i, j]
            value = scale === nothing ? new[i, j] : scale * new[i, j]
            row_change += abs(value - k)
            row_kept += abs(k)
        end
        row_change > row_change_limit(cutoff, row_kept, leak) && return true
    end
    return false
end

# Allowed change of a row with kept norm `norm`: `cutoff` relative to max(1, norm),
# and for the error filter at most A norm^(1-m) (see filter_leak_limit).
@inline function row_change_limit(cutoff, norm, leak)
    limit = cutoff * max(one(norm), norm)
    leak === nothing && return limit
    A, m = leak
    (m > 0 && norm > one(norm)) || return limit
    return min(limit, A * norm^(1 - m))
end

# Estimate local error from the filtered embedded defect for Taylor-Padé,
# from the filtered highest Taylor term for μ-Taylor, or from divided differences.
# The filter (I + (w_1/Q) d J)^Q uses J = ∂f/∂u to damp stiff components
# without the higher Taylor-coefficient Jacobians. For Q = 1 it is M.
function implicit_taylor_error!(integrator, cache, t, dt, uprev, u)
    (; μ, tmp, tmp_embedded, atmp, uintermediate, rhs_embedded, utilde, uprev2, uprev3) = cache
    alg = unwrap_alg(integrator, true)
    (; internalnorm, abstol, reltol) = integrator.opts
    P, Q = _unwrap_val(alg.order), _unwrap_val(alg.order_q)
    if uses_embedded(alg) || (uses_divided3(alg) && integrator.success_iter < 2)
        # Reuse the embedded residual already filtered by the convergence check.
        if cache.embedded_filtered
            copyto!(tmp, tmp_embedded)
        else
            @. tmp = utilde + tmp_embedded - rhs_embedded
            w1 = pade_increment_coefficients(Val(P), Val(Q))[2][2]
            filter_error!(integrator, cache, tmp, Q, w1, dt, Q > 1 ? cache.Mtmp : nothing)
        end
    elseif uses_last_term(alg)
        # tmp_embedded = c_P d^P, d = -h (μ = 1) or (1 - μ)h; filter at -μh.
        w1 = taylor_increment_coefficients(Val(P))[2]
        if μ == one(μ)
            tmp .= tmp_embedded
            filter_error!(integrator, cache, tmp, P, w1, -μ * dt, cache.Mtmp)
        else
            factor = (-μ / (1 - μ))^P
            @. tmp = factor * tmp_embedded
            filter_error!(integrator, cache, tmp, P, w1, -μ * dt, cache.Mtmp)
            @. tmp = tmp - tmp_embedded
        end
    elseif uses_divided3(alg)
        # C h³ y''' with y'''/6 from third divided differences (factor 7/12, as in `Trapezoid`).
        tprev, tprev2 = integrator.tprev, cache.tprev2
        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        dt3 = (t - tprev) * (t - tprev2)
        dt4 = (tprev - tprev2) * (t - tprev2)
        dt5 = t + dt - tprev2
        r = (7 // 12) * 6 * pade_error_constant(P, Q) * dt^3
        @. tmp = r * (
            ((u - uprev) / dt1 - (uprev - uprev2) / dt2) -
                ((uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4)
        ) / dt5
    else
        # dt²/2 y'' from second divided differences (factor 7/12).
        if integrator.success_iter > 0
            tprev = integrator.tprev
            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            r = (7 // 12) * dt^2
            @. tmp = r * ((u - uprev) / dt1 - (uprev - uprev2) / dt2)
        else
            OrdinaryDiffEqCore.set_EEst!(integrator, 1)
            return nothing
        end
    end
    @. atmp = abs(tmp)
    calculate_residuals!(atmp, atmp, uprev, u, abstol, reltol, internalnorm, t)
    OrdinaryDiffEqCore.set_EEst!(integrator, internalnorm(atmp, t))
    return nothing
end

# r ← (I + (w1/degree) d J)^{-degree} r, optionally also for `companion`.
# J is Jf (endpoint ∂c_1/∂u) if given, else the stored ∂c_1/∂u.
function filter_error!(
        integrator, cache, r, degree, w1, distance, Jf = nothing, companion = nothing;
        count = degree
    )
    (; linsolve, linsolve_filter, W, JW, Jk) = cache
    stats = integrator.stats
    if degree == 1
        # Degree one: the filter is M, already factorized.
        dw = linsolve.u
        if cache.direct_lu
            copyto!(dw, _vec(r))
            dense_solve!(cache.J, cache.pivJ, dw)
        else
            linres = checked_dolinsolve(integrator, linsolve; b = _vec(r), linu = dw)
            cache.detached_solution && copyto!(dw, linres.u)
        end
        stats.nsolve += 1
        r .= shaped_like(dw, r)
        return nothing
    end
    alg = unwrap_alg(integrator, true)
    γ = (w1 / degree) * distance
    n = size(W, 1)
    dw = linsolve_filter.u
    right = prepare_filter_rhs!(linsolve_filter.b, r, companion)
    nrhs = companion === nothing ? 1 : 2
    if filter_outdated(cache, γ, filter_cutoff(alg), Jf, filter_leak_limit(alg, cache.κ))
        @inbounds for j in 1:n, i in 1:n
            JW[i, j] = γ * (Jf === nothing ? Jk[i, j, 1] : Jf[i, j])
        end
        matrix_row_norms!(cache.JWrow_norm, JW)
        copyto!(W, JW)
        @inbounds for i in 1:n
            W[i, i] += one(eltype(W))
        end
        cache.filter_valid = true
        stats.nw += 1
        if cache.direct_lu
            dense_lu!(W, cache.pivW)
            copyto!(dw, right)
            dense_solve!(W, cache.pivW, dw)
        else
            linres = checked_dolinsolve(integrator, linsolve_filter; A = W, b = right, linu = dw)
            cache.detached_filter && copyto!(dw, linres.u)
        end
    elseif cache.direct_lu
        copyto!(dw, right)
        dense_solve!(W, cache.pivW, dw)
    else
        linres = checked_dolinsolve(integrator, linsolve_filter; b = right, linu = dw)
        cache.detached_filter && copyto!(dw, linres.u)
    end
    stats.nsolve += nrhs
    copy_filter_result!(r, companion, dw)
    apply_filter_factors!(integrator, cache, r, count - 1, companion)
    return nothing
end

# Apply the kept filter factor `count` times.
@inline function apply_filter_factors!(integrator, cache, r, count, companion = nothing)
    linsolve_filter = cache.linsolve_filter
    dw = linsolve_filter.u
    nrhs = companion === nothing ? 1 : 2
    for _ in 1:count
        right = prepare_filter_rhs!(linsolve_filter.b, r, companion)
        if cache.direct_lu
            copyto!(dw, right)
            dense_solve!(cache.W, cache.pivW, dw)
        else
            linres = checked_dolinsolve(integrator, linsolve_filter; b = right, linu = dw)
            cache.detached_filter && copyto!(dw, linres.u)
        end
        integrator.stats.nsolve += nrhs
        copy_filter_result!(r, companion, dw)
    end
    return nothing
end

@inline prepare_filter_rhs!(b::AbstractVector, r, companion) = _vec(r)
@inline function prepare_filter_rhs!(b::AbstractMatrix, r, companion)
    copyto!(@view(b[:, 1]), _vec(r))
    if companion === nothing
        fill!(@view(b[:, 2]), zero(eltype(b)))
    else
        copyto!(@view(b[:, 2]), _vec(companion))
    end
    return b
end
@inline function copy_filter_result!(r, companion, result::AbstractVector)
    copyto!(_vec(r), result)
    return nothing
end
@inline function copy_filter_result!(r, companion, result::AbstractMatrix)
    copyto!(_vec(r), @view(result[:, 1]))
    companion === nothing || copyto!(_vec(companion), @view(result[:, 2]))
    return nothing
end

# Whether the kept filter differs from γ ∂c_1/∂u by more than `cutoff` in some row
# (see newton_matrix_outdated).
function filter_outdated(cache, γ, cutoff, Jf = nothing, leak = nothing)
    cache.filter_valid || return true
    (; JW, Jk) = cache
    return rows_changed(
        Jf === nothing ? @view(Jk[:, :, 1]) : Jf, JW, cutoff,
        cache.JWrow_norm, cache.row_change, γ, leak
    )
end

newton_state(cache::ImplicitTaylorCache) = cache
newton_state(cache::ImplicitTaylorConstantCache) = cache.stage_cache

# A recovered step still restricts growth, but recovery alone is not a measured
# slow contraction. Record the final attempt's evidence separately for diagnostics.
function record_newton_convergence!(cache, measured, η, alg)
    cache.ηold = η
    cache.newton_rate_measured = measured
    cache.newton_slow = measured && η > fast_convergence_cutoff(alg)
    cache.newton_limit = cache.recovery_since_accept || cache.newton_slow
    return nothing
end

# Scale the previous leak multiplier for growth in |dt|, then raise it to
# the 4/5 power. Without a measured rate, use the current correction’s leak
# directly (Hairer & Wanner, IV.8).
function leak_stop_rate(cache, h)
    η = cache.ηleak
    η < one(η) || return one(η)
    hη = cache.ηleak_h
    scaled = hη > zero(hη) ? η * max(one(h), abs(h) / hη) : η
    return min(max(scaled, eps(η)), one(η))^(4 // 5)
end

# Cap |dt| at 4/5 of the failed size and back off if a cap release just failed.
function record_newton_failure!(state, dt)
    dt = abs(dt)
    state.ηleak = one(state.ηleak)
    state.cap_probe = state.cap_raised ? min(2 * state.cap_probe, 4) : 1
    state.cap_raised = false
    state.cap_streak = 0
    state.dt_cap = dt * convert(typeof(dt), 4 // 5)
    return nothing
end

# Use the controller’s effective factor to select the stage mode. Then limit
# growth after a failed, slow, or recovered solve (q ≥ 1). Release the cap
# after enough fast steps.
function newton_feedback!(state, q, dt, controller_q = q)
    state.stage_auto && stage_mode_feedback!(state, controller_q)
    limited = state.newton_failed || state.newton_limit || state.recovery_since_accept
    if limited
        q = max(q, one(q))
        state.cap_streak = 0
    end
    dt = abs(dt)
    if iszero(state.dt_cap)
        # A fast step after a release confirms it: probe every step again.
        if !limited && state.cap_raised
            state.cap_probe = 1
            state.cap_raised = false
        end
    elseif !limited && dt >= state.dt_cap * convert(typeof(dt), 9 // 10)
        state.cap_streak += 1
        if state.cap_streak >= state.cap_probe
            state.dt_cap = zero(state.dt_cap)
            state.cap_streak = 0
            state.cap_raised = true
        end
    end
    state.newton_failed = false
    state.recovery_since_accept = false
    return q
end

function apply_newton_cap(state, dtnew)
    if !iszero(state.dt_cap) && abs(dtnew) > state.dt_cap
        dtnew = copysign(state.dt_cap, dtnew)
    end
    return dtnew
end

for C in (:PIControllerCache, :PredictiveControllerCache)
    @eval function OrdinaryDiffEqCore.step_accept_controller!(
            integrator, cache::OrdinaryDiffEqCore.$C, alg::ImplicitTaylor, q
        )
        state = newton_state(integrator.cache)
        # Update the error controller’s accepted-step history once.
        dtnew = invoke(
            OrdinaryDiffEqCore.step_accept_controller!,
            Tuple{Any, OrdinaryDiffEqCore.$C, Any, Any}, integrator, cache, alg, q
        )
        controller_q = state.stage_auto ? integrator.dt / dtnew : q
        limited_q = newton_feedback!(state, q, integrator.dt, controller_q)
        if limited_q > q && abs(dtnew) > abs(integrator.dt)
            dtnew = integrator.dt
        end
        return apply_newton_cap(state, dtnew)
    end
end

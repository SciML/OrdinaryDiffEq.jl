function save_idxsinitialize(
        integrator, cache::OrdinaryDiffEqCache,
        ::Type{uType}
    ) where {uType}
    error("This algorithm does not have an initialization function")
end

function loopheader!(integrator)
    # Apply right after iterators / callbacks

    # Accept or reject the step
    if integrator.iter > 0
        if (integrator.opts.adaptive && !integrator.accept_step) ||
                integrator.force_stepfail
            @SciMLMessage(
                lazy"Step rejected: t = $(integrator.t), EEst = $(integrator.EEst)",
                integrator.opts.verbose, :step_rejected
            )
            if integrator.isout
                integrator.dt = integrator.dt * integrator.opts.qmin
            elseif !integrator.force_stepfail
                step_reject_controller!(integrator, integrator.alg)
            end
        else
            @SciMLMessage(
                lazy"Step accepted: t = $(integrator.t), dt = $(integrator.dt), EEst = $(integrator.EEst)",
                integrator.opts.verbose, :step_accepted
            )
            integrator.success_iter += 1
            apply_step!(integrator)
        end
    elseif integrator.u_modified # && integrator.iter == 0
        if integrator.isdae
            DiffEqBase.initialize_dae!(integrator)
        end
        update_uprev!(integrator)
        update_fsal!(integrator)
    end

    integrator.iter += 1
    choose_algorithm!(integrator, integrator.cache)
    fix_dt_at_bounds!(integrator)
    modify_dt_for_tstops!(integrator)
    integrator.force_stepfail = false
    return nothing
end

function apply_step!(integrator)
    update_uprev!(integrator)

    #Update dt if adaptive or if fixed and the dt is allowed to change
    if integrator.opts.adaptive || integrator.dtchangeable
        integrator.dt = integrator.dtpropose
    elseif integrator.dt != integrator.dtpropose && !integrator.dtchangeable
        error("The current setup does not allow for changing dt.")
    end

    update_fsal!(integrator)
    return nothing
end

function update_fsal!(integrator)
    return if has_discontinuity(integrator) &&
            first_discontinuity(integrator) == integrator.tdir * integrator.t
        handle_discontinuities!(integrator)
        get_current_isfsal(integrator.alg, integrator.cache) && reset_fsal!(integrator)
    elseif all_fsal(integrator.alg, integrator.cache) ||
            get_current_isfsal(integrator.alg, integrator.cache)
        if integrator.reeval_fsal || integrator.u_modified ||
                (isdp8(integrator.alg) && !integrator.opts.calck) ||
                (
                only_diagonal_mass_matrix(integrator.alg) &&
                    !integrator.opts.adaptive
            )
            reset_fsal!(integrator)
        else # Do not reeval_fsal, instead copyto! over
            if isinplace(integrator.sol.prob)
                recursivecopy!(integrator.fsalfirst, integrator.fsallast)
            else
                integrator.fsalfirst = integrator.fsallast
            end
        end
    end
end

function last_step_failed(integrator::ODEIntegrator)
    return integrator.last_stepfail && !integrator.opts.adaptive
end

function modify_dt_for_tstops!(integrator)
    return if has_tstop(integrator)
        tdir_t = integrator.tdir * integrator.t
        tdir_tstop = first_tstop(integrator)
        if integrator.opts.adaptive
            integrator.dt = integrator.tdir *
                min(abs(integrator.dt), abs(tdir_tstop - tdir_t)) # step! to the end
        elseif iszero(integrator.dtcache) && integrator.dtchangeable
            integrator.dt = integrator.tdir * abs(tdir_tstop - tdir_t)
        elseif integrator.dtchangeable && !integrator.force_stepfail
            # always try to step! with dtcache, but lower if a tstop
            # however, if force_stepfail then don't set to dtcache, and no tstop worry
            integrator.dt = integrator.tdir *
                min(abs(integrator.dtcache), abs(tdir_tstop - tdir_t)) # step! to the end
        end
    end
end

# Want to extend savevalues! for DDEIntegrator
function savevalues!(integrator::ODEIntegrator, force_save = false, reduce_size = true)
    return _savevalues!(integrator, force_save, reduce_size)
end

function _savevalues!(integrator, force_save, reduce_size)::Tuple{Bool, Bool}
    saved, savedexactly = false, false
    !integrator.opts.save_on && return saved, savedexactly
    tdir_t = integrator.tdir * integrator.t
    saveat = integrator.opts.saveat
    while !isempty(saveat) && first(saveat) <= tdir_t # Perform saveat
        integrator.saveiter += 1
        saved = true
        curt = integrator.tdir * pop!(saveat)
        if curt != integrator.t # If <t, interpolate
            SciMLBase.addsteps!(integrator)
            Θ = (curt - integrator.tprev) / integrator.dt
            val = ode_interpolant(Θ, integrator, integrator.opts.save_idxs, Val{0}) # out of place, but no force copy later
            copyat_or_push!(integrator.sol.t, integrator.saveiter, curt)
            save_val = val
            copyat_or_push!(integrator.sol.u, integrator.saveiter, save_val, false)
            if integrator.alg isa OrdinaryDiffEqCompositeAlgorithm
                copyat_or_push!(
                    integrator.sol.alg_choice, integrator.saveiter,
                    integrator.cache.current
                )
            end
        else # ==t, just save
            if curt == integrator.sol.prob.tspan[2] && !integrator.opts.save_end
                integrator.saveiter -= 1
                continue
            end
            savedexactly = true
            copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
            if integrator.opts.save_idxs === nothing
                copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
            else
                copyat_or_push!(
                    integrator.sol.u, integrator.saveiter,
                    integrator.u[integrator.opts.save_idxs], false
                )
            end
            if isdiscretealg(integrator.alg) || integrator.opts.dense
                integrator.saveiter_dense += 1
                if integrator.opts.dense
                    if integrator.opts.save_idxs === nothing
                        copyat_or_push!(
                            integrator.sol.k, integrator.saveiter_dense,
                            integrator.k
                        )
                    else
                        copyat_or_push!(
                            integrator.sol.k, integrator.saveiter_dense,
                            [k[integrator.opts.save_idxs] for k in integrator.k],
                            false
                        )
                    end
                end
            end
            if integrator.alg isa OrdinaryDiffEqCompositeAlgorithm
                copyat_or_push!(
                    integrator.sol.alg_choice, integrator.saveiter,
                    integrator.cache.current
                )
            end
        end
    end
    if force_save || (
            integrator.opts.save_everystep &&
                (
                isempty(integrator.sol.t) ||
                    (integrator.t !== integrator.sol.t[end]) &&
                    (integrator.opts.save_end || integrator.t !== integrator.sol.prob.tspan[2])
            )
        )
        integrator.saveiter += 1
        saved, savedexactly = true, true
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(
                integrator.sol.u, integrator.saveiter,
                integrator.u[integrator.opts.save_idxs], false
            )
        end
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if isdiscretealg(integrator.alg) || integrator.opts.dense
            integrator.saveiter_dense += 1
            if integrator.opts.dense
                if integrator.opts.save_idxs === nothing
                    copyat_or_push!(
                        integrator.sol.k, integrator.saveiter_dense,
                        integrator.k
                    )
                else
                    copyat_or_push!(
                        integrator.sol.k, integrator.saveiter_dense,
                        [k[integrator.opts.save_idxs] for k in integrator.k],
                        false
                    )
                end
            end
        end
        if integrator.alg isa OrdinaryDiffEqCompositeAlgorithm
            copyat_or_push!(
                integrator.sol.alg_choice, integrator.saveiter,
                integrator.cache.current
            )
        end
    end
    reduce_size && resize!(integrator.k, integrator.kshortsize)
    return saved, savedexactly
end

# Want to extend postamble! for DDEIntegrator
postamble!(integrator::ODEIntegrator) = _postamble!(integrator)

function _postamble!(integrator)
    DiffEqBase.finalize!(integrator.opts.callback, integrator.u, integrator.t, integrator)
    solution_endpoint_match_cur_integrator!(integrator)
    resize!(integrator.sol.t, integrator.saveiter)
    resize!(integrator.sol.u, integrator.saveiter)
    if !(integrator.sol isa DAESolution)
        resize!(integrator.sol.k, integrator.saveiter_dense)
    end
    if integrator.opts.progress
        final_progress(integrator)
    end
    return nothing
end

function final_progress(integrator)
    return @logmsg(
        LogLevel(-1),
        integrator.opts.progress_name,
        _id = integrator.opts.progress_id,
        message = integrator.opts.progress_message(
            integrator.dt, integrator.u,
            integrator.p, integrator.t
        ),
        progress = "done"
    )
end

function solution_endpoint_match_cur_integrator!(integrator)
    return if integrator.opts.save_end &&
            (
            integrator.saveiter == 0 ||
                integrator.sol.t[integrator.saveiter] != integrator.t &&
                (
                (integrator.opts.save_end_user isa Bool && integrator.opts.save_end_user) ||
                    integrator.t ∈ integrator.opts.saveat_cache ||
                    integrator.t == integrator.sol.prob.tspan[2] ||
                    isempty(integrator.opts.saveat_cache)
            )
        )
        integrator.saveiter += 1
        copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
        if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)
        else
            copyat_or_push!(
                integrator.sol.u, integrator.saveiter,
                integrator.u[integrator.opts.save_idxs], false
            )
        end
        if isdiscretealg(integrator.alg) || integrator.opts.dense
            integrator.saveiter_dense += 1
            if integrator.opts.dense
                if integrator.opts.save_idxs === nothing
                    copyat_or_push!(
                        integrator.sol.k, integrator.saveiter_dense,
                        integrator.k
                    )
                else
                    copyat_or_push!(
                        integrator.sol.k, integrator.saveiter_dense,
                        [k[integrator.opts.save_idxs] for k in integrator.k],
                        false
                    )
                end
            end
        end
        if integrator.alg isa OrdinaryDiffEqCompositeAlgorithm
            copyat_or_push!(
                integrator.sol.alg_choice, integrator.saveiter,
                integrator.cache.current
            )
        end
        SciMLBase.save_final_discretes!(integrator, integrator.opts.callback)
    end
end

# Want to extend loopfooter! for DDEIntegrator
loopfooter!(integrator::ODEIntegrator) = _loopfooter!(integrator)

function _loopfooter!(integrator)
    # Carry-over from callback
    # This is set to true if u_modified requires callback FSAL reset
    # But not set to false when reset so algorithms can check if reset occurred
    integrator.reeval_fsal = false
    integrator.u_modified = false
    integrator.do_error_check = true
    ttmp = integrator.t + integrator.dt
    if integrator.force_stepfail
        if integrator.opts.adaptive
            post_newton_controller!(integrator, integrator.alg)
        elseif integrator.last_stepfail
            return
        end
        integrator.last_stepfail = true
        integrator.accept_step = false
    elseif integrator.opts.adaptive
        q = stepsize_controller!(integrator, integrator.alg)
        integrator.isout = integrator.opts.isoutofdomain(integrator.u, integrator.p, ttmp)
        integrator.accept_step = (
            !integrator.isout &&
                accept_step_controller(
                integrator,
                integrator.opts.controller
            )
        ) ||
            (
            integrator.opts.force_dtmin &&
                abs(integrator.dt) <= timedepentdtmin(integrator)
        )
        if integrator.accept_step # Accept
            increment_accept!(integrator.stats)
            integrator.last_stepfail = false
            dtnew = DiffEqBase.value(
                step_accept_controller!(
                    integrator,
                    integrator.alg,
                    q
                )
            ) *
                oneunit(integrator.dt)
            integrator.tprev = integrator.t
            integrator.t = fixed_t_for_floatingpoint_error!(integrator, ttmp)
            calc_dt_propose!(integrator, dtnew)
            handle_callbacks!(integrator)
        else # Reject
            increment_reject!(integrator.stats)
        end
    elseif !integrator.opts.adaptive #Not adaptive
        increment_accept!(integrator.stats)
        integrator.tprev = integrator.t
        integrator.t = fixed_t_for_floatingpoint_error!(integrator, ttmp)
        integrator.last_stepfail = false
        integrator.accept_step = true
        integrator.dtpropose = integrator.dt
        handle_callbacks!(integrator)
    end
    if integrator.opts.progress && integrator.iter % integrator.opts.progress_steps == 0
        log_step!(
            integrator.opts.progress_name, integrator.opts.progress_id,
            integrator.opts.progress_message, integrator.dt, integrator.u,
            integrator.p, integrator.t, integrator.sol.prob.tspan
        )
    end

    # Take value because if t is dual then maxeig can be dual
    if integrator.cache isa CompositeCache
        cur_eigen_est = integrator.opts.internalnorm(
            DiffEqBase.value(integrator.eigen_est),
            integrator.t
        )
        cur_eigen_est > integrator.stats.maxeig &&
            (integrator.stats.maxeig = cur_eigen_est)
    end
    return nothing
end

function increment_accept!(stats)
    return stats.naccept += 1
end

function increment_reject!(stats)
    return stats.nreject += 1
end

function log_step!(progress_name, progress_id, progress_message, dt, u, p, t, tspan)
    t1, t2 = tspan
    return @logmsg(
        LogLevel(-1), progress_name,
        _id = progress_id,
        message = progress_message(dt, u, p, t),
        progress = (t - t1) / (t2 - t1)
    )
end

function fixed_t_for_floatingpoint_error!(integrator, ttmp)
    return if has_tstop(integrator)
        tstop = integrator.tdir * first_tstop(integrator)
        if abs(ttmp - tstop) <
                100eps(float(max(integrator.t, tstop) / oneunit(integrator.t))) *
                oneunit(integrator.t)
            tstop
        else
            ttmp
        end
    else
        ttmp
    end
end

# Use a generated function to call apply_callback! in a type-stable way
@generated function apply_ith_callback!(
        integrator,
        time, upcrossing, event_idx, cb_idx,
        callbacks::NTuple{
            N,
            Union{
                ContinuousCallback,
                VectorContinuousCallback,
            },
        }
    ) where {N}
    ex = quote
        throw(BoundsError(callbacks, cb_idx))
    end
    for i in 1:N
        # N.B: doing this as an explicit if (return) else (rest of expression)
        # means that LLVM compiles this into a switch.
        # This seemingly isn't the case with just if (return) end (rest of expression)
        ex = quote
            if (cb_idx == $i)
                return DiffEqBase.apply_callback!(
                    integrator, callbacks[$i], time,
                    upcrossing, event_idx
                )
            else
                $ex
            end
        end
    end
    return ex
end

function handle_callbacks!(integrator)
    discrete_callbacks = integrator.opts.callback.discrete_callbacks
    continuous_callbacks = integrator.opts.callback.continuous_callbacks
    atleast_one_callback = false

    continuous_modified = false
    discrete_modified = false
    saved_in_cb = false
    if !(continuous_callbacks isa Tuple{})
        time, upcrossing,
            event_occurred,
            event_idx,
            idx,
            counter = DiffEqBase.find_first_continuous_callback(
            integrator,
            continuous_callbacks...
        )
        if event_occurred
            integrator.event_last_time = idx
            integrator.vector_event_last_time = event_idx
            continuous_modified,
                saved_in_cb = apply_ith_callback!(
                integrator,
                time, upcrossing,
                event_idx,
                idx,
                continuous_callbacks
            )
        else
            integrator.event_last_time = 0
            integrator.vector_event_last_time = 1
        end
    end
    if !integrator.force_stepfail && !(discrete_callbacks isa Tuple{})
        discrete_modified,
            saved_in_cb = DiffEqBase.apply_discrete_callback!(
            integrator,
            discrete_callbacks...
        )
    end
    if !saved_in_cb
        savevalues!(integrator)
    end

    integrator.u_modified = continuous_modified | discrete_modified
    integrator.reeval_fsal && handle_callback_modifiers!(integrator) # Hook for DDEs to add discontinuities
    return nothing
end

function update_uprev!(integrator)
    if alg_extrapolates(integrator.alg)
        if isinplace(integrator.sol.prob)
            recursivecopy!(integrator.uprev2, integrator.uprev)
        else
            integrator.uprev2 = integrator.uprev
        end
    end
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.uprev, integrator.u)
        if integrator.alg isa DAEAlgorithm
            recursivecopy!(integrator.duprev, integrator.du)
        end
    else
        integrator.uprev = integrator.u
        if integrator.alg isa DAEAlgorithm
            integrator.duprev = integrator.du
        end
    end
    return nothing
end

handle_discontinuities!(integrator) = pop_discontinuity!(integrator)

function calc_dt_propose!(integrator, dtnew)
    dtnew = if has_dtnew_modification(integrator.alg) &&
            integrator.opts.adaptive && (integrator.iter >= 1)
        dtnew_modification(integrator, integrator.alg, dtnew)
    else
        dtnew
    end
    dtpropose = integrator.tdir * min(abs(integrator.opts.dtmax), abs(dtnew))
    dtpropose = integrator.tdir * max(abs(dtpropose), timedepentdtmin(integrator))
    integrator.dtpropose = dtpropose
    return nothing
end

function fix_dt_at_bounds!(integrator)
    if integrator.tdir > 0
        integrator.dt = min(integrator.opts.dtmax, integrator.dt)
    else
        integrator.dt = max(integrator.opts.dtmax, integrator.dt)
    end
    dtmin = timedepentdtmin(integrator)
    if integrator.tdir > 0
        integrator.dt = max(integrator.dt, dtmin)
    else
        integrator.dt = min(integrator.dt, dtmin)
    end
    return nothing
end

function handle_tstop!(integrator)
    if has_tstop(integrator)
        tdir_t = integrator.tdir * integrator.t
        tdir_tstop = first_tstop(integrator)
        if tdir_t == tdir_tstop
            while tdir_t == tdir_tstop #remove all redundant copies
                res = pop_tstop!(integrator)
                has_tstop(integrator) ? (tdir_tstop = first_tstop(integrator)) : break
            end
            integrator.just_hit_tstop = true
        elseif tdir_t > tdir_tstop
            if !integrator.dtchangeable
                SciMLBase.change_t_via_interpolation!(
                    integrator,
                    integrator.tdir *
                        pop_tstop!(integrator), Val{true}
                )
                integrator.just_hit_tstop = true
            else
                error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
            end
        end
    end
    return nothing
end

handle_callback_modifiers!(integrator::ODEIntegrator) = nothing

function reset_fsal!(integrator)
    # Under these conditions, these algorithms are not FSAL anymore
    increment_nf!(integrator.stats, 1)

    # Ignore DAEs but they already re-ran initialization
    # Mass matrix DAEs do need to reset FSAL if available
    return if !(integrator.sol.prob isa DAEProblem)
        if ismutablecache(integrator.cache)
            integrator.f(integrator.fsalfirst, integrator.u, integrator.p, integrator.t)
        else
            integrator.fsalfirst = integrator.f(integrator.u, integrator.p, integrator.t)
        end
    end

    # Do not set false here so it can be checked in the algorithm
    # integrator.reeval_fsal = false
end

function nlsolve_f(f, alg::OrdinaryDiffEqAlgorithm)
    return f isa SplitFunction && issplit(alg) ? f.f1 : f
end
nlsolve_f(f, alg::DAEAlgorithm) = f
function nlsolve_f(integrator::ODEIntegrator)
    return nlsolve_f(integrator.f, unwrap_alg(integrator, true))
end

function (integrator::ODEIntegrator)(
        t, ::Type{deriv} = Val{0};
        idxs = nothing
    ) where {deriv}
    return current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::ODEIntegrator)(
        val::AbstractArray, t::Union{Number, AbstractArray},
        ::Type{deriv} = Val{0}; idxs = nothing
    ) where {deriv}
    return current_interpolant!(val, t, integrator, idxs, deriv)
end

has_discontinuity(integrator) = !isempty(integrator.opts.d_discontinuities)
first_discontinuity(integrator) = first(integrator.opts.d_discontinuities)
pop_discontinuity!(integrator) = pop!(integrator.opts.d_discontinuities)

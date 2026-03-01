# Noise interface functions — no-ops when W/P are nothing (pure ODE).
# StochasticDiffEq extends these with methods for NoiseProcess types.
accept_noise!(::Nothing, args...) = nothing
reject_noise!(::Nothing, args...) = nothing
save_noise!(::Nothing) = nothing
noise_curt(::Nothing) = nothing
is_noise_saveable(::Nothing) = false

# Trait: does the solution type support dense output k-array storage?
# True for ODESolution (has sol.k), false for RODESolution/DAESolution (no sol.k).
@inline _has_ks(integrator) = hasfield(typeof(integrator.sol), :k)

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
        if (!integrator.force_stepfail) &&
                (
                !integrator.opts.adaptive || integrator.accept_step ||
                    isaposteriori(integrator.alg)
            )
            # ACCEPT
            @SciMLMessage(
                lazy"Step accepted: t = $(integrator.t), dt = $(integrator.dt), EEst = $(integrator.EEst)",
                integrator.opts.verbose, :step_accepted
            )
            integrator.success_iter += 1
            apply_step!(integrator)
        elseif (
                integrator.opts.adaptive && !integrator.accept_step &&
                    !isaposteriori(integrator.alg)
            ) ||
                integrator.force_stepfail
            # REJECT
            handle_step_rejection!(integrator)
        end
    elseif integrator.u_modified # && integrator.iter == 0
        on_u_modified_at_init!(integrator)
    end

    integrator.iter += 1
    choose_algorithm!(integrator, integrator.cache)
    fix_dt_at_bounds!(integrator)
    modify_dt_for_tstops!(integrator)
    integrator.force_stepfail = false
    return nothing
end

# Handles step rejection in loopheader: adjust dt, reject noise, and call post_step_reject!.
function handle_step_rejection!(integrator)
    @SciMLMessage(
        lazy"Step rejected: t = $(integrator.t), EEst = $(integrator.EEst)",
        integrator.opts.verbose, :step_rejected
    )
    if integrator.isout
        integrator.dt = integrator.dt * integrator.opts.qmin
    elseif !integrator.force_stepfail
        step_reject_controller!(integrator, integrator.alg)
    end
    # Noise rejection (no-op when W/P are nothing for pure ODE)
    if !isnothing(integrator.W)
        fix_dt_at_bounds!(integrator)
        modify_dt_for_tstops!(integrator)
        reject_noise!(integrator.W, integrator.dt, integrator.u, integrator.p)
        reject_noise!(integrator.P, integrator.dt, integrator.u, integrator.p)
        integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))
    end
    return post_step_reject!(integrator)
end

# Called after step rejection handling. Override for DDE discontinuity handling.
post_step_reject!(integrator) = nothing

# Called at iter==0 when u was modified by callbacks during init.
# For SDE: isdae=false skips DAE re-init; isfsal=false makes update_fsal! a no-op.
function on_u_modified_at_init!(integrator)
    if integrator.isdae
        DiffEqBase.initialize_dae!(integrator)
    end
    update_uprev!(integrator)
    return update_fsal!(integrator)
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

    # Shorten dt to hit the next tstop after update_fsal!, which for DDEs calls
    # handle_discontinuities! using integrator.dt to track propagated discontinuities
    # in the interval [t, t+dt]. Must come after update_fsal! but before
    # noise acceptance so that SDE noise sees the tstop-adjusted dt.
    modify_dt_for_tstops!(integrator)

    # Noise acceptance (no-op when W/P are nothing for pure ODE)
    accept_noise!(integrator.W, integrator.dt, integrator.u, integrator.p, true)
    accept_noise!(integrator.P, integrator.dt, integrator.u, integrator.p, true)
    if !isnothing(integrator.W)
        integrator.dt = integrator.W.dt  # RSWM readback
        integrator.sqdt = @fastmath integrator.tdir * sqrt(abs(integrator.dt))
    end

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

# Accessor functions for tstop flag fields with fallbacks for non-ODE integrators
# (e.g. DDEIntegrator in DelayDiffEq.jl which doesn't have these fields)
_get_next_step_tstop(integrator::ODEIntegrator) = integrator.next_step_tstop
_get_next_step_tstop(integrator) = false

function _set_tstop_flag!(integrator::ODEIntegrator, is_tstop::Bool, target = nothing)
    integrator.next_step_tstop = is_tstop
    if is_tstop && target !== nothing
        integrator.tstop_target = target
    end
    return nothing
end
_set_tstop_flag!(integrator, is_tstop::Bool, target = nothing) = nothing

_get_tstop_target(integrator::ODEIntegrator) = integrator.tstop_target

function modify_dt_for_tstops!(integrator)
    return if has_tstop(integrator)
        tdir_t = integrator.tdir * integrator.t
        tdir_tstop = first_tstop(integrator)
        distance_to_tstop = abs(tdir_tstop - tdir_t)

        if integrator.opts.adaptive
            original_dt = abs(integrator.dt)
            integrator.dtpropose = original_dt
            if original_dt < distance_to_tstop
                _set_tstop_flag!(integrator, false)
            else
                _set_tstop_flag!(
                    integrator, true, integrator.tdir * tdir_tstop)
            end
            integrator.dt = integrator.tdir * min(original_dt, distance_to_tstop)
        elseif iszero(integrator.dtcache) && integrator.dtchangeable
            integrator.dt = integrator.tdir * distance_to_tstop
            _set_tstop_flag!(
                integrator, true, integrator.tdir * tdir_tstop)
        elseif integrator.dtchangeable && !integrator.force_stepfail
            # always try to step! with dtcache, but lower if a tstop
            # however, if force_stepfail then don't set to dtcache, and no tstop worry
            if abs(integrator.dtcache) < distance_to_tstop
                _set_tstop_flag!(integrator, false)
            else
                _set_tstop_flag!(
                    integrator, true, integrator.tdir * tdir_tstop)
            end
            integrator.dt = integrator.tdir *
                min(abs(integrator.dtcache), distance_to_tstop)
        else
            _set_tstop_flag!(integrator, false)
        end
    else
        _set_tstop_flag!(integrator, false)
    end
end

function handle_tstop_step!(integrator)
    return if integrator.t isa AbstractFloat && abs(integrator.dt) < eps(abs(integrator.t))
        # Skip perform_step! entirely for tiny dt
        integrator.accept_step = true
    else
        perform_step!(integrator, integrator.cache)
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
            Θ = (curt - integrator.tprev) / integrator.dt
            val = interp_at_saveat(Θ, integrator, integrator.opts.save_idxs, Val{0})
            copyat_or_push!(integrator.sol.t, integrator.saveiter, curt)
            save_val = val
            copyat_or_push!(integrator.sol.u, integrator.saveiter, save_val, false)
            if is_composite_algorithm(integrator.alg)
                copyat_or_push!(
                    integrator.sol.alg_choice, integrator.saveiter,
                    integrator.cache.current
                )
            end
        else # ==t, just save
            if skip_saveat_at_tspan_end(integrator, curt)
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
            save_dense_at_t!(integrator)
            if is_composite_algorithm(integrator.alg)
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
                    (integrator.t !== integrator.sol.t[end] || iszero(integrator.dt)) &&
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
        save_dense_at_t!(integrator)
        if is_composite_algorithm(integrator.alg)
            copyat_or_push!(
                integrator.sol.alg_choice, integrator.saveiter,
                integrator.cache.current
            )
        end
    end
    post_savevalues!(integrator, reduce_size)
    return saved, savedexactly
end

# Interpolation at saveat points.
# ODE (opts.dense = true): polynomial interpolation via addsteps!/ode_interpolant.
# SDE (opts.dense = false): linear interpolation between uprev and u.
function interp_at_saveat(Θ, integrator, idxs, deriv)
    if integrator.opts.dense
        SciMLBase.addsteps!(integrator)
        return ode_interpolant(Θ, integrator, idxs, deriv)
    else
        return linear_interpolant(Θ, integrator, idxs, deriv)
    end
end

# Linear interpolation: (1 - Θ) * uprev + Θ * u
@inline function linear_interpolant(Θ, integrator, idxs::Nothing, ::Type{Val{0}})
    return @. (1 - Θ) * integrator.uprev + Θ * integrator.u
end
@inline function linear_interpolant(Θ, integrator, idxs, ::Type{Val{0}})
    return @. (1 - Θ) * integrator.uprev[idxs] + Θ * integrator.u[idxs]
end
@inline function linear_interpolant(Θ, integrator, idxs::Nothing, ::Type{Val{1}})
    return @. (integrator.u - integrator.uprev) / integrator.dt
end
@inline function linear_interpolant(Θ, integrator, idxs, ::Type{Val{1}})
    return @. (integrator.u[idxs] - integrator.uprev[idxs]) / integrator.dt
end

# Skip saveat at tspan end when save_end=false.
# ODE: skip saving at tspan[2] when save_end=false.
# SDE: always save at explicit saveat times (never skip).
function skip_saveat_at_tspan_end(integrator, curt)
    return isnothing(integrator.W) && curt == integrator.sol.prob.tspan[2] &&
        !integrator.opts.save_end
end

# Save dense output when saving at exact time t.
# Only stores k-array data when the solution supports it (_has_ks).
function save_dense_at_t!(integrator)
    return if (isdiscretealg(integrator.alg) || integrator.opts.dense) && _has_ks(integrator)
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
end

# Cleanup after savevalues: resize k for dense output storage.
# No-op when solution lacks k-array storage (SDE/RODE).
function post_savevalues!(integrator, reduce_size)
    return reduce_size && integrator.opts.dense && _has_ks(integrator) &&
        resize!(integrator.k, integrator.kshortsize)
end

# Want to extend postamble! for DDEIntegrator
postamble!(integrator::ODEIntegrator) = _postamble!(integrator)

function _postamble!(integrator)
    DiffEqBase.finalize!(integrator.opts.callback, integrator.u, integrator.t, integrator)
    solution_endpoint_match_cur_integrator!(integrator)
    resize!(integrator.sol.t, integrator.saveiter)
    resize!(integrator.sol.u, integrator.saveiter)
    finalize_solution_storage!(integrator)
    if integrator.opts.progress
        final_progress(integrator)
    end
    return nothing
end

# Finalize solution storage in postamble: resize arrays, save noise.
function finalize_solution_storage!(integrator)
    sizehint!(integrator.sol.t, integrator.saveiter)
    sizehint!(integrator.sol.u, integrator.saveiter)
    # Dense output arrays (only when solution has k-array storage)
    if integrator.opts.dense && _has_ks(integrator) && !(integrator.sol isa DAESolution)
        resize!(integrator.sol.k, integrator.saveiter_dense)
        sizehint!(integrator.sol.k, integrator.saveiter_dense)
    end
    # Noise finalization (SDE only)
    if !isnothing(integrator.W) && noise_curt(integrator.W) != integrator.t
        accept_noise!(integrator.W, integrator.dt, integrator.u, integrator.p, false)
        accept_noise!(integrator.P, integrator.dt, integrator.u, integrator.p, false)
    end
    if is_noise_saveable(integrator.W) && !integrator.W.save_everystep
        save_noise!(integrator.W)
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
        if (isdiscretealg(integrator.alg) || integrator.opts.dense) && _has_ks(integrator)
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
        if is_composite_algorithm(integrator.alg)
            copyat_or_push!(
                integrator.sol.alg_choice, integrator.saveiter,
                integrator.cache.current
            )
        end
        finalize_endpoint!(integrator)
    end
end

# Called at the end of solution_endpoint_match_cur_integrator!: save final discretes.
function finalize_endpoint!(integrator)
    return SciMLBase.save_final_discretes!(integrator, integrator.opts.callback)
end

# Want to extend loopfooter! for DDEIntegrator
loopfooter!(integrator::ODEIntegrator) = _loopfooter!(integrator)

function _loopfooter!(integrator)
    loopfooter_reset!(integrator)
    integrator.do_error_check = true
    ttmp = integrator.t + integrator.dt
    if integrator.force_stepfail
        if integrator.opts.adaptive
            handle_force_stepfail!(integrator)
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
            integrator.tprev = integrator.t

            if _get_next_step_tstop(integrator)
                # Step controller dt is overly pessimistic, since dt = time to tstop.
                # Restore the original dt so the controller proposes a reasonable next step.
                integrator.dt = integrator.dtpropose
            end
            integrator.t = fixed_t_for_tstop_error!(integrator, ttmp)

            dtnew = DiffEqBase.value(
                step_accept_controller!(
                    integrator,
                    integrator.alg,
                    q
                )
            ) *
                oneunit(integrator.dt)
            calc_dt_propose!(integrator, dtnew)
            handle_callbacks!(integrator)
        else # Reject
            increment_reject!(integrator.stats)
        end
    elseif !integrator.opts.adaptive #Not adaptive
        increment_accept!(integrator.stats)
        integrator.tprev = integrator.t
        integrator.t = fixed_t_for_tstop_error!(integrator, ttmp)
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
    if is_composite_cache(integrator.cache)
        cur_eigen_est = integrator.opts.internalnorm(
            DiffEqBase.value(integrator.eigen_est),
            integrator.t
        )
        cur_eigen_est > integrator.stats.maxeig &&
            (integrator.stats.maxeig = cur_eigen_est)
    end
    return nothing
end

# Trait: is this a composite algorithm cache? Override to include SDE composite caches.
is_composite_cache(cache) = cache isa CompositeCache

# Trait: is this a composite algorithm? Override to include SDE composite algorithms.
is_composite_algorithm(alg) = alg isa OrdinaryDiffEqCompositeAlgorithm

# Reset integrator flags at the start of loopfooter.
# For SDE, reeval_fsal is always false so resetting is a no-op.
function loopfooter_reset!(integrator)
    # Carry-over from callback
    # This is set to true if u_modified requires callback FSAL reset
    # But not set to false when reset so algorithms can check if reset occurred
    integrator.reeval_fsal = false
    return integrator.u_modified = false
end

# Handle force_stepfail in adaptive mode: reduce dt after Newton failure.
# post_newton_controller! does dt = dt / failfactor, which works for both ODE and SDE.
handle_force_stepfail!(integrator) = post_newton_controller!(integrator, integrator.alg)

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

function fixed_t_for_tstop_error!(integrator, ttmp)
    if _get_next_step_tstop(integrator)
        _set_tstop_flag!(integrator, false)
        return _get_tstop_target(integrator)
    else
        return ttmp
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
    on_callbacks_complete!(integrator)
    return nothing
end

# Called after all callbacks have been applied.
# ODE/FSAL: trigger FSAL re-evaluation and DDE discontinuity handling.
# SDE: Poisson rate update when u is modified.
function on_callbacks_complete!(integrator)
    if isfsal(integrator.alg)
        integrator.reeval_fsal && handle_callback_modifiers!(integrator)
    elseif integrator.u_modified && !isnothing(integrator.W)
        integrator.do_error_check = false
        handle_callback_modifiers!(integrator)
    end
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

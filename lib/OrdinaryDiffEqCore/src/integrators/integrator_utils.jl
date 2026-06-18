# Noise interface functions — no-ops when W/P are nothing (pure ODE).
# StochasticDiffEq extends these with methods for NoiseProcess types.
accept_noise!(::Nothing, args...) = nothing
reject_noise!(::Nothing, args...) = nothing
save_noise!(::Nothing) = nothing
noise_curt(::Nothing) = nothing
is_noise_saveable(::Nothing) = false
reinit_noise!(::Nothing, dt) = nothing

# Noise field accessors — safe for any integrator type.
# ODEIntegrator has W/P/sqdt; other integrators (DDEIntegrator) don't.
@inline _get_W(integrator) = hasfield(typeof(integrator), :W) ? getfield(integrator, :W) : nothing
@inline _get_P(integrator) = hasfield(typeof(integrator), :P) ? getfield(integrator, :P) : nothing

# Trait: does the integrator+solution support dense output k-array storage?
# True for ODEIntegrator (has integrator.k and sol.k), false for SDEIntegrator
# (no integrator.k) and RODESolution/DAESolution (no sol.k).
@inline _has_ks(integrator) = hasfield(typeof(integrator), :k) && hasfield(typeof(integrator.sol), :k)

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
                lazy"Step accepted: t = $(integrator.t), dt = $(integrator.dt), EEst = $(get_EEst(integrator))",
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
    elseif integrator.derivative_discontinuity # && integrator.iter == 0
        on_derivative_discontinuity_at_init!(integrator)
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
        lazy"Step rejected: t = $(integrator.t), EEst = $(get_EEst(integrator))",
        integrator.opts.verbose, :step_rejected
    )
    if integrator.isout
        integrator.dt = integrator.dt * get_qmin(integrator)
    elseif !integrator.force_stepfail
        step_reject_controller!(integrator, integrator.alg)
    end
    # Noise rejection (no-op when W/P are nothing for pure ODE)
    W = _get_W(integrator)
    if !isnothing(W)
        fix_dt_at_bounds!(integrator)
        modify_dt_for_tstops!(integrator)
        reject_noise!(W, integrator.dt, integrator.u, integrator.p)
        reject_noise!(_get_P(integrator), integrator.dt, integrator.u, integrator.p)
        integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))
    end
    return post_step_reject!(integrator)
end

# Called after step rejection handling. Override for DDE discontinuity handling.
post_step_reject!(integrator) = nothing

# Called at iter==0 when u was modified by callbacks during init.
# For SDE: isdae=false skips DAE re-init; isfsal=false makes update_fsal! a no-op.
function on_derivative_discontinuity_at_init!(integrator)
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
    W = _get_W(integrator)
    accept_noise!(W, integrator.dt, integrator.u, integrator.p, true)
    accept_noise!(_get_P(integrator), integrator.dt, integrator.u, integrator.p, true)
    if !isnothing(W)
        integrator.dt = W.dt  # RSWM readback
        integrator.sqdt = @fastmath integrator.tdir * sqrt(abs(integrator.dt))
    end

    return nothing
end

function update_fsal!(integrator)
    if has_discontinuity(integrator) &&
            first_discontinuity(integrator) == integrator.tdir * integrator.t
        handle_discontinuities!(integrator)
        shift_past_discontinuity!(integrator)
        get_current_isfsal(integrator.alg, integrator.cache) && reset_fsal!(integrator)
    elseif all_fsal(integrator.alg, integrator.cache) ||
            get_current_isfsal(integrator.alg, integrator.cache)
        if integrator.reeval_fsal || integrator.derivative_discontinuity ||
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
    return nothing
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
    if has_tstop(integrator)
        tdir_t = integrator.tdir * integrator.t
        tdir_tstop = first_tstop(integrator)
        distance_to_tstop = abs(tdir_tstop - tdir_t)
        # Floating-point tolerance so that a dt whose nominal value matches
        # distance_to_tstop to within rounding still triggers the tstop
        # branch.  Without this, accumulated `t + dt + dt + …` can drift
        # just past the last tstop and produce a spurious micro-step.
        tstop_tol = if integrator.t isa AbstractFloat && isfinite(tdir_tstop) &&
                isfinite(integrator.t)
            100 * eps(
                float(
                    max(abs(integrator.t), abs(tdir_tstop)) /
                        oneunit(integrator.t)
                )
            ) * oneunit(integrator.t)
        else
            zero(distance_to_tstop)
        end

        if integrator.opts.adaptive
            original_dt = abs(integrator.dt)
            integrator.dtpropose = integrator.tdir * original_dt
            if original_dt + tstop_tol < distance_to_tstop
                _set_tstop_flag!(integrator, false)
            else
                _set_tstop_flag!(
                    integrator, true, integrator.tdir * tdir_tstop
                )
            end
            integrator.dt = integrator.tdir * min(original_dt, distance_to_tstop)
        elseif iszero(integrator.dtcache) && integrator.dtchangeable
            integrator.dt = integrator.tdir * distance_to_tstop
            _set_tstop_flag!(
                integrator, true, integrator.tdir * tdir_tstop
            )
        elseif integrator.dtchangeable && !integrator.force_stepfail
            # always try to step! with dtcache, but lower if a tstop
            # however, if force_stepfail then don't set to dtcache, and no tstop worry
            if abs(integrator.dtcache) + tstop_tol < distance_to_tstop
                _set_tstop_flag!(integrator, false)
            else
                _set_tstop_flag!(
                    integrator, true, integrator.tdir * tdir_tstop
                )
            end
            integrator.dt = integrator.tdir *
                min(abs(integrator.dtcache), distance_to_tstop)
        else
            _set_tstop_flag!(integrator, false)
        end
    else
        _set_tstop_flag!(integrator, false)
    end
    return nothing
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
# ODE: polynomial interpolation via addsteps!/ode_interpolant (always available,
#   regardless of opts.dense which only controls post-solve k-array storage).
# SDE: linear interpolation between uprev and u.
function interp_at_saveat(Θ, integrator, idxs, ::Type{deriv}) where {deriv}
    if isnothing(_get_W(integrator))
        # ODE/DDE: polynomial interpolation
        SciMLBase.addsteps!(integrator)
        return ode_interpolant(Θ, integrator, idxs, deriv)
    else
        # SDE: linear interpolation
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
    return isnothing(_get_W(integrator)) && curt == integrator.sol.prob.tspan[2] &&
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
    return reduce_size && _has_ks(integrator) &&
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
    W = _get_W(integrator)
    if !isnothing(W) && noise_curt(W) != integrator.t
        accept_noise!(W, integrator.dt, integrator.u, integrator.p, false)
        accept_noise!(_get_P(integrator), integrator.dt, integrator.u, integrator.p, false)
    end
    if is_noise_saveable(W) && !W.save_everystep
        save_noise!(W)
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
                integrator.alg
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
    # This is set to true if derivative_discontinuity requires callback FSAL reset
    # But not set to false when reset so algorithms can check if reset occurred
    integrator.reeval_fsal = false
    return integrator.derivative_discontinuity = false
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

# overrides this with a method that calls calc_J to get a fresh Jacobian.
get_fresh_jacobian(integrator, cache) = cache.J

function SciMLBase.log_instability(integrator::ODEIntegrator)
    W = _get_W(integrator)
    u = integrator.u
    u0 = integrator.sol.prob.u0

    # state analysis: NaN/Inf components, and components that have blown up
    nan_inf_idxs = findall(!isfinite, u)
    blown_idxs = Int[]
    if length(u) == length(u0)
        for i in eachindex(u)
            ref = max(abs(u0[i]), oneunit(eltype(u)))
            abs(u[i]) > 1.0e6 * ref && push!(blown_idxs, i)
        end
    end

    # jacobian analysis over rows and columns for large values
    jac = if W !== nothing && hasproperty(W, :J)
        #rosenbrock
        W.J
    elseif hasproperty(integrator.cache, :J)
        #radau
        get_fresh_jacobian(integrator, integrator.cache)
    elseif hasproperty(integrator.cache, :nlsolver) &&
            hasproperty(integrator.cache.nlsolver.cache, :J)
        #BDF
        integrator.cache.nlsolver.cache.J
    else
        nothing
    end

    bad_entries = nothing
    singularity_rows = nothing
    singularity_cols = nothing
    if jac !== nothing
        rows = Set{Int}()
        cols = Set{Int}()
        entries = Tuple{Int, Int, eltype(jac)}[]
        _find_large_jac_entries!(rows, cols, entries, jac)

        # keep only entries within 10 orders of magnitude of the largest finite entry,
        # plus any non-finite entries. filters out large-but-normal model parameters 
        max_finite = 0.0
        for (_, _, v) in entries
            if isfinite(v)
                max_finite = max(max_finite, abs(v))
            end
        end
        cutoff = max_finite * 1e-10
        filter!(t -> !isfinite(t[3]) || abs(t[3]) >= cutoff, entries) #only keep those vals within 1e10 of max or inf/nan
        sort!(entries, by = t -> (!isfinite(t[3]), abs(t[3])), rev = true)

        # derive rows and columns from remaining entries
        row_set = Set{Int}()
        col_set = Set{Int}()
        for (i, j, _) in entries
            push!(row_set, i)
            push!(col_set, j)
        end
        bad_entries = entries
        singularity_rows = sort!(collect(row_set))
        singularity_cols = sort!(collect(col_set))
    end

    # trace diagnostics to symbolic system if present
    f = integrator.sol.prob.f
    sys = (hasproperty(f, :sys) && f.sys !== nothing) ? f.sys : nothing
    sym_eqs = (sys !== nothing && hasfield(typeof(sys), :eqs)) ? getfield(sys, :eqs) : nothing
    sym_vars = (sys !== nothing && hasfield(typeof(sys), :unknowns)) ? getfield(sys, :unknowns) : nothing

    # diagnostic message construction
    diagnostic = String[]
    if !isempty(nan_inf_idxs) #state vars
        if u isa AbstractArray
            n_nan = length(nan_inf_idxs)
            n_total = length(u)
            if n_nan == n_total
                push!(diagnostic, "All $n_total state variables are non-finite (NaN/Inf)")
            elseif n_nan > 3
                push!(diagnostic, "$n_nan of $n_total state variables are non-finite (NaN/Inf): indices $nan_inf_idxs")
            else
                for i in nan_inf_idxs
                    push!(diagnostic, "u[$i] = $(u[i]) is non-finite (NaN/Inf)")
                end
            end
        else
            push!(diagnostic, "u = $u is non-finite (NaN/Inf), suggesting blow-up or a NaN in the RHS")
        end
    elseif !isempty(blown_idxs)
        if u isa AbstractArray
            for i in blown_idxs
                push!(diagnostic, "u[$i] = $(@sprintf("%.4g", u[i])) has grown >1e6× its initial value")
            end
        else
            push!(diagnostic, "u = $(@sprintf("%.4g", u)) has grown >1e6× its initial value")
        end
    end

    if bad_entries !== nothing && !isempty(bad_entries) #Jacobian analysis
        has_nonfinite = false
        has_large = false
        for (_, _, v) in bad_entries
            isfinite(v) ? (has_large = true) : (has_nonfinite = true)
        end
        entry_desc = if has_nonfinite && has_large
            "non-finite and large"
        elseif has_nonfinite
            "non-finite"
        else
            "unusually large"
        end

        example_strs = String[]
        for (i, j, v) in first(bad_entries, 5)
            push!(example_strs, "J[$i,$j] = $(@sprintf("%.4g", v))")
        end
        push!(diagnostic, "Jacobian row(s) $singularity_rows have $entry_desc entries (e.g. $(join(example_strs, ", "))), suggesting a singularity in those equation(s)")
        if sym_eqs !== nothing
            for row in singularity_rows
                if row <= length(sym_eqs)
                    push!(diagnostic, "  row $row corresponds to equation: $(sym_eqs[row])") #trace rows back to symbolic eqs
                end
            end
        end
        # jac cols
        if !isempty(singularity_cols)
            push!(diagnostic, "Jacobian column(s) $singularity_cols have $entry_desc entries, suggesting those state component(s) are diverging")
            if sym_vars !== nothing
                for col in singularity_cols
                    if col <= length(sym_vars)
                        push!(diagnostic, "  col $col corresponds to variable: $(sym_vars[col])") #trace cols back to symbolic vars
                    end
                end
            end
        end
    end
    diagnostic = isempty(diagnostic) ? "" : "\nDiagnostics:\n" * join(diagnostic, "\n") * "."
    return diagnostic
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

    integrator.derivative_discontinuity = continuous_modified | discrete_modified
    on_callbacks_complete!(integrator)
    return nothing
end

# Called after all callbacks have been applied.
# ODE/FSAL: trigger FSAL re-evaluation and DDE discontinuity handling.
# SDE: Poisson rate update when u is modified.
function on_callbacks_complete!(integrator)
    if isfsal(integrator.alg)
        integrator.reeval_fsal && handle_callback_modifiers!(integrator)
    elseif integrator.derivative_discontinuity && !isnothing(_get_W(integrator))
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

"""
    shift_past_discontinuity!(integrator)

Advance `integrator.t` by one ULP in the `integrator.tdir` direction. Called
right after `handle_discontinuities!` so that subsequent RHS evaluations —
including the FSAL re-evaluation that follows in `update_fsal!` — take place
on the post-discontinuity side of `t`-dependent branches in the user's `f`.

By convention, a `d_discontinuities` entry at `t_d` marks the vector field as
right-discontinuous there: `f` evaluated at `t_d` is the "old" regime, and
`f` at `nextfloat(t_d)` is the "new" regime. This pairs with user code
written as `if t > t_d; new; else; old; end` and lets starting-time
discontinuities (`t_d == t0`) work by advancing forward into the tspan.

The state `u` is continuous across `t_d` (only the vector field changes), so
it is left untouched. For non-`AbstractFloat` time types (e.g. `Rational`),
there is no ULP to shift to and the call is a no-op.
"""
@inline function shift_past_discontinuity!(integrator)
    _shift_past_discontinuity!(integrator, integrator.t, integrator.tdir)
    return nothing
end
@inline function _shift_past_discontinuity!(integrator, t::AbstractFloat, tdir)
    integrator.t = tdir > 0 ? nextfloat(t) : prevfloat(t)
    return nothing
end
@inline _shift_past_discontinuity!(integrator, t, tdir) = nothing

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
function nlsolve_f(f, alg::StochasticDiffEqAlgorithm)
    return f isa SplitSDEFunction && issplit(alg) ? f.f1 : f
end
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

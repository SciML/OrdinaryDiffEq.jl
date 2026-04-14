function SciMLBase.__solve(
        prob::Union{SciMLBase.AbstractDDEProblem, AbstractSDDEProblem},
        alg::AbstractMethodOfStepsAlgorithm, args...;
        kwargs...
    )
    integrator = SciMLBase.__init(prob, alg, args...; kwargs...)
    DiffEqBase.solve!(integrator)
    return integrator.sol
end

# Compile-time detection for OrdinaryDiffEqCore version compatibility
# Detects if DEOptions has typeof(verbose) as a type parameter (PR #2895)
"""
    _count_deoptions_typeparams()

Count the number of type parameters in OrdinaryDiffEqCore.DEOptions.
Used for compile-time detection of API changes.
"""
function _count_deoptions_typeparams()
    count = 0
    T = OrdinaryDiffEqCore.DEOptions
    while T isa UnionAll
        count += 1
        T = T.body
    end
    return count
end

"""
    DEOPTIONS_HAS_VERBOSE_TYPEPARAM

Compile-time constant that is `true` if OrdinaryDiffEqCore.DEOptions includes
`typeof(verbose)` as a type parameter (OrdinaryDiffEqCore >= 1.37.0 or with PR #2895).
Old version has 20 parameters, new version has 21 parameters (adds typeof(verbose)).
"""
const DEOPTIONS_HAS_VERBOSE_TYPEPARAM = _count_deoptions_typeparams() >= 21

function SciMLBase.__init(
        prob::Union{SciMLBase.AbstractDDEProblem, AbstractSDDEProblem},
        alg::AbstractMethodOfStepsAlgorithm,
        timeseries_init = (),
        ts_init = (),
        ks_init = ();
        saveat = (),
        tstops = (),
        d_discontinuities = (),
        save_idxs = nothing,
        save_everystep = isempty(saveat),
        save_on = true,
        save_start = save_everystep || isempty(saveat) ||
            saveat isa Number || prob.tspan[1] in saveat,
        save_end = nothing,
        save_discretes = true,
        save_noise = save_everystep,
        callback = nothing,
        dense = save_everystep && isempty(saveat),
        calck = (callback !== nothing && callback != CallbackSet()) || # Empty callback
            dense, # and no dense output
        dt = zero(eltype(prob.tspan)),
        dtmin = DiffEqBase.prob2dtmin(prob; use_end_time = false),
        dtmax = eltype(prob.tspan)(prob.tspan[end] - prob.tspan[1]),
        force_dtmin = false,
        adaptive = isadaptive(alg),
        gamma = OrdinaryDiffEqCore.gamma_default(alg.alg),
        abstol = nothing,
        reltol = nothing,
        qmin = OrdinaryDiffEqCore.qmin_default(alg.alg),
        qmax = OrdinaryDiffEqCore.qmax_default(alg.alg),
        qsteady_min = OrdinaryDiffEqCore.qsteady_min_default(alg.alg),
        qsteady_max = OrdinaryDiffEqCore.qsteady_max_default(alg.alg),
        qoldinit = isadaptive(alg) ? 1 // 10^4 : 0,
        controller = nothing,
        fullnormalize = true,
        failfactor = 2,
        beta1 = nothing,
        beta2 = nothing,
        maxiters = adaptive ? 1000000 : typemax(Int),
        internalnorm = DiffEqBase.ODE_DEFAULT_NORM,
        internalopnorm = opnorm,
        isoutofdomain = DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN,
        unstable_check = DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
        verbose = DEVerbosity(),
        timeseries_errors = true,
        dense_errors = false,
        advance_to_tstop = false,
        stop_at_next_tstop = false,
        initialize_save = true,
        progress = false,
        progress_steps = 1000,
        progress_name = "DDE",
        progress_message = DiffEqBase.ODE_DEFAULT_PROG_MESSAGE,
        progress_id = :DelayDiffEq,
        userdata = nothing,
        allow_extrapolation = OrdinaryDiffEqCore.alg_extrapolates(alg),
        initialize_integrator = true,
        alias_u0 = false,
        seed = UInt64(0),
        # keyword arguments for DDEs
        discontinuity_interp_points::Int = 10,
        discontinuity_abstol = eltype(prob.tspan)(1 // Int64(10)^12),
        discontinuity_reltol = 0,
        initializealg = DDEDefaultInit(),
        kwargs...
    )
    is_stochastic = prob isa AbstractSDDEProblem

    if haskey(kwargs, :initial_order)
        @warn "initial_order has been deprecated. Please specify order_discontinuity_t0 in the DDEProblem/SDDEProblem instead."
        order_discontinuity_t0::Int = kwargs[:initial_order]
    else
        order_discontinuity_t0 = is_stochastic ? Int(prob.order_discontinuity_t0) : prob.order_discontinuity_t0
    end

    # Handle verbose argument: convert Bool or AbstractVerbosityPreset to DEVerbosity
    if verbose isa Bool
        if verbose
            verbose_spec = DEVerbosity()
        else
            verbose_spec = DEVerbosity(None())
        end
    elseif verbose isa AbstractVerbosityPreset
        verbose_spec = DEVerbosity(verbose)
    else
        verbose_spec = verbose
    end

    if !is_stochastic && alg.alg isa CompositeAlgorithm && alg.alg.choice_function isa AutoSwitch
        auto = alg.alg.choice_function
        alg = MethodOfSteps(
            CompositeAlgorithm(
                alg.alg.algs,
                OrdinaryDiffEqCore.AutoSwitchCache(
                    0, 0,
                    auto.nonstiffalg,
                    auto.stiffalg,
                    auto.stiffalgfirst,
                    auto.maxstiffstep,
                    auto.maxnonstiffstep,
                    auto.nonstifftol,
                    auto.stifftol,
                    auto.dtfac,
                    auto.stiffalgfirst,
                    auto.switch_max
                )
            )
        )
    end

    if haskey(kwargs, :minimal_solution)
        @warn "minimal_solution is ignored"
    end

    if !isempty(saveat) && dense
        @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
    end

    progress && @logmsg(-1, progress_name, _id = progress_id, progress = 0)

    isdae = prob.f.mass_matrix !== I && !(prob.f.mass_matrix isa Tuple) &&
        ArrayInterface.issingular(prob.f.mass_matrix)

    # unpack problem
    (; f, u0, h, tspan, p, neutral, constant_lags, dependent_lags) = prob

    # determine type and direction of time
    tType = eltype(tspan)
    t0 = first(tspan)
    tdir = sign(last(tspan) - t0)
    tTypeNoUnits = typeof(one(tType))

    # Allow positive dtmax, but auto-convert
    dtmax > zero(dtmax) && tdir < zero(tdir) && (dtmax *= tdir)

    # no fixed-point iterations for constrained algorithms,
    # and thus `dtmax` should match minimal lag
    if isconstrained(alg) && has_constant_lags(prob)
        min_lag = minimum(abs, constant_lags)
        old_dtmax = abs(dtmax)
        dtmax = tdir * min(old_dtmax, min_lag)
        if min_lag < old_dtmax
            @SciMLMessage(
                lazy"Constrained algorithm: limiting dtmax from $old_dtmax to $min_lag (minimum lag)",
                verbose_spec, :constrained_step
            )
        end
    end

    # get absolute and relative tolerances
    abstol_internal = get_abstol(u0, tspan, alg.alg; abstol = abstol)
    reltol_internal = get_reltol(u0, tspan, alg.alg; reltol = reltol)

    # get rate prototype
    rate_prototype = rate_prototype_of(u0, tspan)

    # compute noise_rate_prototype for SDDE
    if is_stochastic
        _noise_rate_prototype = prob.noise_rate_prototype
        if is_diagonal_noise(prob)
            noise_rate_prototype = rate_prototype
        elseif _noise_rate_prototype !== nothing
            noise_rate_prototype = copy(_noise_rate_prototype)
        else
            noise_rate_prototype = nothing
        end
    else
        noise_rate_prototype = nothing
    end

    # get states (possibly different from the ODE integrator!)
    u, uprev,
        uprev2 = u_uprev_uprev2(
        u0, alg;
        alias_u0 = alias_u0,
        adaptive = adaptive,
        allow_extrapolation = allow_extrapolation,
        calck = calck
    )
    uEltypeNoUnits = recursive_unitless_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    # get the differential vs algebraic variables
    differential_vars = if is_stochastic
        OrdinaryDiffEqCore.get_differential_vars(f, u)
    elseif prob isa DAEProblem
        prob.differential_vars
    else
        OrdinaryDiffEqCore.get_differential_vars(f, u)
    end

    # create a history function
    history = build_history_function(
        prob, alg, rate_prototype, reltol_internal,
        differential_vars;
        dt = dt, dtmin = dtmin, calck = false,
        adaptive = adaptive, internalnorm = internalnorm
    )
    f_with_history = if is_stochastic
        SDEFunctionWrapper(f, history)
    else
        ODEFunctionWrapper(f, history)
    end

    # ── Noise creation (SDDE only) ─────────────────────────────────────
    if is_stochastic
        W, P, sqdt = _create_sdde_noise(
            prob, alg.alg, u0, t0, tType(dt), tdir, noise_rate_prototype,
            save_noise, seed, isinplace(prob), isadaptive(alg)
        )
    else
        W = nothing
        P = nothing
        sqdt = nothing
    end

    # initialize output arrays of the solution
    save_idxs,
        saved_subsystem = SciMLBase.get_save_idxs_and_saved_subsystem(prob, save_idxs)

    k = typeof(rate_prototype)[]
    ts, timeseries,
        ks = solution_arrays(
        u, tspan, rate_prototype;
        timeseries_init = timeseries_init,
        ts_init = ts_init,
        ks_init = ks_init,
        save_idxs = save_idxs,
        save_start = save_start,
        is_stochastic = is_stochastic
    )

    # build cache
    ode_integrator = history.integrator
    cache = if is_stochastic
        dW = W.dW
        dZ = W.dZ
        _sde_alg_cache(
            alg.alg, prob, u, dW, dZ, p,
            rate_prototype, noise_rate_prototype,
            nothing, uEltypeNoUnits, uBottomEltypeNoUnits,
            tTypeNoUnits, uprev, f_with_history, t0,
            tType(dt), Val{isinplace(prob)},
            DiffEqBase.DEVerbosity()
        )
    else
        OrdinaryDiffEqCore.alg_cache(
            alg.alg, u, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2,
            f_with_history, t0, zero(tType), reltol_internal, p,
            calck,
            Val(isinplace(prob)), DiffEqBase.DEVerbosity()
        )
    end

    # separate statistics of the integrator and the history
    stats = SciMLBase.DEStats(0)

    # create solution
    alg_choice = iscomposite(alg) ? Int[] : nothing
    id = OrdinaryDiffEqCore.InterpolationData(
        f_with_history, timeseries, ts, ks,
        alg_choice, dense, cache, differential_vars, false
    )
    sol = if is_stochastic
        SciMLBase.build_solution(
            prob, alg.alg, ts, timeseries;
            dense = dense, k = ks, interp = id, saved_subsystem = saved_subsystem,
            alg_choice = id.alg_choice, calculate_error = false,
            stats = stats, W = W
        )
    else
        SciMLBase.build_solution(
            prob, alg.alg, ts, timeseries;
            dense = dense, k = ks, interp = id, saved_subsystem = saved_subsystem,
            alg_choice = id.alg_choice, calculate_error = false,
            stats = stats
        )
    end

    # retrieve time stops, time points at which solutions is saved, and discontinuities
    tstops_internal = OrdinaryDiffEqCore.initialize_tstops(
        tType, tstops, d_discontinuities,
        tspan
    )
    saveat_internal = OrdinaryDiffEqCore.initialize_saveat(tType, saveat, tspan)
    d_discontinuities_internal = OrdinaryDiffEqCore.initialize_d_discontinuities(
        Discontinuity{
            tType,
            Int,
        },
        d_discontinuities,
        tspan
    )

    maximum_order = OrdinaryDiffEqCore.alg_maximum_order(alg)
    tstops_propagated,
        d_discontinuities_propagated = initialize_tstops_d_discontinuities_propagated(
        tType,
        tstops,
        d_discontinuities,
        tspan,
        order_discontinuity_t0,
        maximum_order,
        constant_lags,
        neutral
    )

    # reserve capacity for the solution
    _sizehint_solution!(
        sol, alg, tspan, tstops_internal, saveat_internal;
        save_everystep = save_everystep, adaptive = adaptive, dt = tType(dt),
        dtmin = dtmin, internalnorm = internalnorm
    )

    # create array of tracked discontinuities
    # used to find propagated discontinuities with callbacks and to keep track of all
    # passed discontinuities
    tracked_discontinuities = Discontinuity{tType, Int}[]
    if order_discontinuity_t0 ≤ maximum_order
        push!(tracked_discontinuities, Discontinuity(tdir * t0, order_discontinuity_t0))
    end

    # Create set of callbacks and its cache
    callback_set, callback_cache = callback_set_and_cache(prob, callback)

    # separate options of integrator and of dummy ODE integrator since ODE integrator always saves
    # every step and every index (necessary for history function)
    QT = tTypeNoUnits <: Integer ? typeof(qmin) : typeof(internalnorm(u, t0))

    # Setting up the step size controller
    if (beta1 !== nothing || beta2 !== nothing) && controller !== nothing
        throw(ArgumentError("Setting both the legacy PID parameters `beta1, beta2 = $((beta1, beta2))` and the `controller = $controller` is not allowed."))
    end

    if (beta1 !== nothing || beta2 !== nothing)
        message = "Providing the legacy PID parameters `beta1, beta2` is deprecated. Use the keyword argument `controller` instead."
        Base.depwarn(message, :init)
        Base.depwarn(message, :solve)
    end

    if controller === nothing
        controller = OrdinaryDiffEqCore.default_controller(
            alg.alg, cache,
            convert(QT, qoldinit)::QT, beta1,
            beta2
        )
    end

    save_end_user = save_end
    save_end = save_end === nothing ?
        save_everystep || isempty(saveat) || saveat isa Number ||
        prob.tspan[2] in saveat : save_end

    # Compute delta for SDE algorithms (used by Milstein-type methods).
    # For pure DDE, delta is nothing (unused).
    delta = is_stochastic ?
        convert(recursive_unitless_bottom_eltype(u), 1 // 1) : nothing

    # Construct DEOptions using full constructor (with delta and save_noise fields)
    opts = OrdinaryDiffEqCore.DEOptions{
        typeof(abstol_internal), typeof(reltol_internal),
        QT, tType, typeof(controller),
        typeof(internalnorm), typeof(internalopnorm),
        typeof(save_end_user),
        typeof(callback_set),
        typeof(isoutofdomain),
        typeof(progress_message), typeof(unstable_check),
        typeof(tstops_internal),
        typeof(d_discontinuities_internal), typeof(userdata),
        typeof(save_idxs),
        typeof(maxiters), typeof(tstops),
        typeof(saveat), typeof(d_discontinuities), typeof(verbose_spec),
        typeof(delta),
    }(
        maxiters,
        save_everystep,
        adaptive,
        abstol_internal,
        reltol_internal,
        QT(gamma),
        QT(qmax),
        QT(qmin),
        QT(qsteady_max),
        QT(qsteady_min),
        QT(qoldinit),
        QT(failfactor),
        tType(dtmax),
        tType(dtmin),
        controller,
        internalnorm,
        internalopnorm,
        save_idxs,
        tstops_internal,
        saveat_internal,
        d_discontinuities_internal,
        tstops,
        saveat,
        d_discontinuities,
        userdata,
        progress,
        progress_steps,
        progress_name,
        progress_message,
        progress_id,
        timeseries_errors,
        dense_errors,
        delta,
        dense,
        save_on,
        save_start,
        save_end,
        save_noise,
        save_discretes,
        save_end_user,
        callback_set,
        isoutofdomain,
        unstable_check,
        verbose_spec,
        calck,
        force_dtmin,
        advance_to_tstop,
        stop_at_next_tstop
    )

    # create fixed point solver
    fpsolver = build_fpsolver(
        alg, alg.fpsolve, u, uEltypeNoUnits, uBottomEltypeNoUnits,
        Val(isinplace(prob))
    )

    # initialize indices of u(t) and u(tprev) in the dense history
    prev_idx = 1
    prev2_idx = 1

    # create integrator combining the new defined problem function with history
    # information, the new solution, the parameters of the ODE integrator, and
    # parameters of fixed-point iteration
    # do not initialize fsalfirst and fsallast
    # rate/state = (state/time)/state = 1/t units, internalnorm drops units
    eigen_est = inv(one(tType))
    tprev = t0
    dtcache = tType(dt)
    dtpropose = tType(dt)
    iter = 0
    kshortsize = 0
    reeval_fsal = false
    u_modified = false
    EEst = QT(1)
    just_hit_tstop = false
    do_error_check = true
    isout = false
    accept_step = false
    force_stepfail = false
    last_stepfail = false
    event_last_time = 0
    vector_event_last_time = 1
    last_event_error = zero(uBottomEltypeNoUnits)
    dtchangeable = OrdinaryDiffEqCore.isdtchangeable(alg.alg)
    q11 = QT(1)
    success_iter = 0
    erracc = QT(1)
    dtacc = tType(1)

    fsalfirst, fsallast = OrdinaryDiffEqCore.get_fsalfirstlast(cache, rate_prototype)
    integrator = DDEIntegrator{
        typeof(alg.alg), isinplace(prob), typeof(u0), tType,
        typeof(p),
        typeof(eigen_est), QT, typeof(tdir), typeof(k), typeof(sol),
        typeof(f_with_history), typeof(cache),
        typeof(ode_integrator), typeof(fpsolver),
        typeof(opts), typeof(discontinuity_abstol),
        typeof(discontinuity_reltol), typeof(history),
        typeof(tstops_propagated),
        typeof(d_discontinuities_propagated),
        typeof(fsalfirst),
        typeof(last_event_error), typeof(callback_cache),
        typeof(differential_vars), typeof(initializealg),
        typeof(W), typeof(P), typeof(sqdt), Nothing,
    }(
        sol, u, k,
        t0,
        tType(dt),
        f_with_history,
        p,
        uprev,
        uprev2,
        tprev,
        prev_idx,
        prev2_idx,
        fpsolver,
        order_discontinuity_t0,
        tracked_discontinuities,
        discontinuity_interp_points,
        discontinuity_abstol,
        discontinuity_reltol,
        tstops_propagated,
        d_discontinuities_propagated,
        alg.alg,
        dtcache,
        dtchangeable,
        dtpropose,
        tdir,
        eigen_est,
        EEst,
        QT(qoldinit),
        q11,
        erracc,
        dtacc,
        success_iter,
        iter,
        length(ts),
        length(ts),
        cache,
        callback_cache,
        kshortsize,
        force_stepfail,
        last_stepfail,
        just_hit_tstop,
        do_error_check,
        event_last_time,
        vector_event_last_time,
        last_event_error,
        accept_step,
        isout,
        reeval_fsal,
        u_modified,
        isdae,
        opts,
        stats,
        history,
        differential_vars,
        ode_integrator, fsalfirst, fsallast, initializealg,
        W, P, sqdt, nothing,
    )

    # initialize DDE integrator
    if initialize_integrator
        SciMLBase.initialize_dae!(integrator)
        initialize_solution!(integrator)
        OrdinaryDiffEqCore.initialize_callbacks!(integrator, initialize_save)
        if save_on && save_start
            SciMLBase.save_discretes_if_enabled!(integrator, opts.callback; skip_duplicates = true)
        end
        DiffEqBase.initialize!(integrator)
    end

    # take care of time step dt = 0 and dt with incorrect sign
    OrdinaryDiffEqCore.handle_dt!(integrator)

    # After handle_dt! may have changed dt (via auto_dt_reset!), update noise state.
    # This mirrors what OrdinaryDiffEqCore's ODE __init does for SDE integrators.
    if !isnothing(integrator.W)
        OrdinaryDiffEqCore.modify_dt_for_tstops!(integrator)
        integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))
        integrator.W.dt = integrator.dt
    end

    # Starting-time discontinuity handling: mirrors OrdinaryDiffEqCore's __init.
    OrdinaryDiffEqCore.handle_starting_time_discontinuity!(integrator)

    return integrator
end

function DiffEqBase.solve!(integrator::DDEIntegrator)
    (; tdir, opts, sol) = integrator
    (; tstops) = opts

    # step over all stopping time points, similar to solving with ODE integrators
    @inbounds while !isempty(tstops)
        while tdir * integrator.t < first(tstops)
            # apply step or adapt step size
            OrdinaryDiffEqCore.loopheader!(integrator)

            # abort integration following same criteria as for ODEs:
            # maxiters exceeded, dt <= dtmin, integration unstable
            SciMLBase.check_error!(integrator) == ReturnCode.Success || return sol

            # calculate next step
            OrdinaryDiffEqCore.perform_step!(integrator)

            # calculate proposed next step size, handle callbacks, and update solution
            OrdinaryDiffEqCore.loopfooter!(integrator)

            isempty(tstops) && break
        end

        # remove hit or passed stopping time points
        OrdinaryDiffEqCore.handle_tstop!(integrator)
    end

    # clean up solution
    SciMLBase.postamble!(integrator)

    f = sol.prob.f

    if SciMLBase.has_analytic(f)
        SciMLBase.calculate_solution_errors!(
            sol;
            timeseries_errors = opts.timeseries_errors,
            dense_errors = opts.dense_errors
        )
    end
    sol.retcode == ReturnCode.Default || return sol

    return integrator.sol = SciMLBase.solution_new_retcode(sol, ReturnCode.Success)
end

function OrdinaryDiffEqCore.initialize_callbacks!(
        integrator::DDEIntegrator,
        initialize_save = true
    )
    callbacks = integrator.opts.callback
    prob = integrator.sol.prob

    # set up additional initial values of newly created DDE integrator
    # (such as fsalfirst) and its callbacks

    integrator.u_modified = true

    u_modified = initialize!(callbacks, integrator.u, integrator.t, integrator)

    # if the user modifies u, we need to fix previous values before initializing
    # FSAL in order for the starting derivatives to be correct
    if u_modified
        if isinplace(prob)
            recursivecopy!(integrator.uprev, integrator.u)
        else
            integrator.uprev = integrator.u
        end

        if OrdinaryDiffEqCore.alg_extrapolates(integrator.alg)
            if isinplace(prob)
                recursivecopy!(integrator.uprev2, integrator.uprev)
            else
                integrator.uprev2 = integrator.uprev
            end
        end

        # update heap of discontinuities
        # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
        add_next_discontinuities!(integrator, 0, integrator.t)

        # reset this as it is now handled so the integrators should proceed as normal
        reeval_internals_due_to_modification!(integrator, Val{false})

        if initialize_save &&
                (
                any((c) -> c.save_positions[2], callbacks.discrete_callbacks) ||
                    any((c) -> c.save_positions[2], callbacks.continuous_callbacks)
            )
            savevalues!(integrator, true)
        end
    end

    # reset this as it is now handled so the integrators should proceed as normal
    return integrator.u_modified = false
end

function initialize_tstops_d_discontinuities_propagated(
        ::Type{T}, tstops,
        d_discontinuities, tspan,
        order_discontinuity_t0,
        alg_maximum_order,
        constant_lags, neutral
    ) where {T}
    # create heaps for propagated discontinuities and corresponding time stops
    tstops_propagated = BinaryMinHeap{T}()
    d_discontinuities_propagated = BinaryMinHeap{Discontinuity{T, Int}}()

    # add discontinuities and time stops propagated from initial discontinuity
    if constant_lags !== nothing && !isempty(constant_lags) &&
            order_discontinuity_t0 ≤ alg_maximum_order
        sizehint!(tstops_propagated, length(constant_lags))
        sizehint!(d_discontinuities_propagated, length(constant_lags))

        t0, tf = tspan
        tdir = sign(tf - t0)
        maxlag = tdir * (tf - t0)
        next_order = neutral ? order_discontinuity_t0 : order_discontinuity_t0 + 1
        for lag in constant_lags
            if tdir * lag < maxlag
                t = tdir * (t0 + lag)
                push!(tstops_propagated, t)

                d = Discontinuity{T, Int}(t, next_order)
                push!(d_discontinuities_propagated, d)
            end
        end
    end

    return tstops_propagated, d_discontinuities_propagated
end

# Override check_prob_alg_pairing: MethodOfSteps wrapping an SDE algorithm can solve SDDEProblem
function DiffEqBase.check_prob_alg_pairing(prob::SDDEProblem, alg::AbstractMethodOfStepsAlgorithm)
    # MethodOfSteps wrapping an SDE algorithm is valid for SDDEProblem
    if !(alg.alg isa SDEAlgUnion)
        throw(SciMLBase.ProblemSolverPairingError(prob, alg))
    end
    if isdefined(prob, :u0) && SciMLBase.eltypedual(prob.u0) &&
            !SciMLBase.isautodifferentiable(alg)
        throw(SciMLBase.DirectAutodiffError())
    end
    return nothing
end

# SDE caches don't have fsalfirst/fsallast (SDE algorithms are never FSAL)
OrdinaryDiffEqCore.get_fsalfirstlast(cache::StochasticDiffEqConstantCache, u) = (zero(u), zero(u))
OrdinaryDiffEqCore.get_fsalfirstlast(cache::StochasticDiffEqMutableCache, u) = (zero(u), zero(u))

## Note: StochasticDiffEqCore already defines initialize!(integrator, cache::StochasticDiffEqCache)

struct DDEDefaultInit <: SciMLBase.DAEInitializationAlgorithm end

function SciMLBase.initialize_dae!(integrator::DDEIntegrator, initializealg = integrator.initializealg)
    return OrdinaryDiffEqCore._initialize_dae!(
        integrator, integrator.sol.prob, initializealg,
        Val(DiffEqBase.isinplace(integrator.sol.prob))
    )
end

function OrdinaryDiffEqCore._initialize_dae!(integrator::DDEIntegrator, prob, ::DDEDefaultInit, isinplace)
    return if SciMLBase.has_initializeprob(prob.f)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.OverrideInit(), isinplace)
    else
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.CheckInit(), isinplace)
    end
end

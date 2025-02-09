function DiffEqBase.__solve(
        prob::Union{DiffEqBase.AbstractODEProblem,
            DiffEqBase.AbstractDAEProblem},
        alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, args...;
        kwargs...)
    integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
    solve!(integrator)
    integrator.sol
end

function DiffEqBase.__init(
        prob::Union{DiffEqBase.AbstractODEProblem,
            DiffEqBase.AbstractDAEProblem},
        alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm},
        timeseries_init = (),
        ts_init = (),
        ks_init = (),
        recompile::Type{Val{recompile_flag}} = Val{true};
        saveat = (),
        tstops = (),
        d_discontinuities = (),
        save_idxs = nothing,
        save_everystep = isempty(saveat),
        save_on = true,
        save_start = save_everystep || isempty(saveat) ||
                         saveat isa Number || prob.tspan[1] in saveat,
        save_end = nothing,
        callback = nothing,
        dense = save_everystep && isempty(saveat) && !default_linear_interpolation(prob, alg),
        calck = (callback !== nothing && callback !== CallbackSet()) ||
                    (dense) || !isempty(saveat), # and no dense output
        dt = isdiscretealg(alg) && isempty(tstops) ?
             eltype(prob.tspan)(1) : eltype(prob.tspan)(0),
        dtmin = eltype(prob.tspan)(0),
        dtmax = eltype(prob.tspan)((prob.tspan[end] - prob.tspan[1])),
        force_dtmin = false,
        adaptive = anyadaptive(alg),
        gamma = gamma_default(alg),
        abstol = nothing,
        reltol = nothing,
        qmin = qmin_default(alg),
        qmax = qmax_default(alg),
        qsteady_min = qsteady_min_default(alg),
        qsteady_max = qsteady_max_default(alg),
        beta1 = nothing,
        beta2 = nothing,
        qoldinit = anyadaptive(alg) ? 1 // 10^4 : 0,
        controller = nothing,
        fullnormalize = true,
        failfactor = 2,
        maxiters = anyadaptive(alg) ? 1000000 : typemax(Int),
        internalnorm = ODE_DEFAULT_NORM,
        internalopnorm = LinearAlgebra.opnorm,
        isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
        unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
        verbose = true,
        timeseries_errors = true,
        dense_errors = false,
        advance_to_tstop = false,
        stop_at_next_tstop = false,
        initialize_save = true,
        progress = false,
        progress_steps = 1000,
        progress_name = "ODE",
        progress_message = ODE_DEFAULT_PROG_MESSAGE,
        progress_id = :OrdinaryDiffEq,
        userdata = nothing,
        allow_extrapolation = alg_extrapolates(alg),
        initialize_integrator = true,
        alias = ODEAliasSpecifier(),
        initializealg = DefaultInit(),
        kwargs...) where {recompile_flag}
    if prob isa DiffEqBase.AbstractDAEProblem && alg isa OrdinaryDiffEqAlgorithm
        error("You cannot use an ODE Algorithm with a DAEProblem")
    end

    if prob isa DiffEqBase.AbstractODEProblem && alg isa DAEAlgorithm
        error("You cannot use an DAE Algorithm with a ODEProblem")
    end

    if prob isa DiffEqBase.ODEProblem
        if !(prob.f isa DiffEqBase.DynamicalODEFunction) && alg isa PartitionedAlgorithm
            error("You can not use a solver designed for partitioned ODE with this problem. Please choose a solver suitable for your problem")
        end
    end

    if prob.f isa DynamicalODEFunction && prob.f.mass_matrix isa Tuple
        if any(mm != I for mm in prob.f.mass_matrix)
            error("This solver is not able to use mass matrices. For compatible solvers see https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/")
        end
    elseif !(prob isa DiscreteProblem) &&
           !(prob isa DiffEqBase.AbstractDAEProblem) &&
           !is_mass_matrix_alg(alg) &&
           prob.f.mass_matrix != I
        error("This solver is not able to use mass matrices. For compatible solvers see https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/")
    end

    if alg isa OrdinaryDiffEqRosenbrockAdaptiveAlgorithm &&
       # https://github.com/SciML/OrdinaryDiffEq.jl/pull/2079 fixes this for Rosenbrock23 and 32
       !only_diagonal_mass_matrix(alg) &&
       prob.f.mass_matrix isa AbstractMatrix &&
       all(isequal(0), prob.f.mass_matrix)
        # technically this should also warn for zero operators but those are hard to check for
        if (dense || !isempty(saveat)) && verbose
            @warn("Rosenbrock methods on equations without differential states do not bound the error on interpolations.")
        end
    end

    if only_diagonal_mass_matrix(alg) &&
       prob.f.mass_matrix isa AbstractMatrix &&
       !isdiag(prob.f.mass_matrix)
        throw(ArgumentError("$(typeof(alg).name.name) only works with diagonal mass matrices. Please choose a solver suitable for your problem (e.g. Rodas5P)"))
    end

    if !isempty(saveat) && dense
        @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
    end

    progress && @logmsg(LogLevel(-1), progress_name, _id=progress_id, progress=0)

    tType = eltype(prob.tspan)
    tspan = prob.tspan
    tdir = sign(tspan[end] - tspan[1])

    t = tspan[1]

    if (((!(alg isa OrdinaryDiffEqAdaptiveAlgorithm) &&
          !(alg isa OrdinaryDiffEqCompositeAlgorithm) &&
          !(alg isa DAEAlgorithm)) || !adaptive || !isadaptive(alg)) &&
        dt == tType(0) && isempty(tstops)) && dt_required(alg)
        throw(ArgumentError("Fixed timestep methods require a choice of dt or choosing the tstops"))
    end
    if !isadaptive(alg) && adaptive
        throw(ArgumentError("Fixed timestep methods can not be run with adaptive=true"))
    end

    isdae = alg isa DAEAlgorithm || (!(prob isa DiscreteProblem) &&
             prob.f.mass_matrix != I &&
             !(prob.f.mass_matrix isa Tuple) &&
             ArrayInterface.issingular(prob.f.mass_matrix))
    if alg isa CompositeAlgorithm && alg.choice_function isa AutoSwitch
        auto = alg.choice_function
        _alg = CompositeAlgorithm(alg.algs,
            AutoSwitchCache(0, 0,
                auto.nonstiffalg,
                auto.stiffalg,
                auto.stiffalgfirst,
                auto.maxstiffstep,
                auto.maxnonstiffstep,
                auto.nonstifftol,
                auto.stifftol,
                auto.dtfac,
                auto.stiffalgfirst,
                auto.switch_max, 0))
    else
        _alg = alg
    end

    use_old_kwargs = haskey(kwargs,:alias_u0) || haskey(kwargs,:alias_du0)

    if use_old_kwargs
        aliases = ODEAliasSpecifier()
        if haskey(kwargs, :alias_u0)
            message = "`alias_u0` keyword argument is deprecated, to set `alias_u0`,
            please use an ODEAliasSpecifier, e.g. `solve(prob, alias = ODEAliasSpecifier(alias_u0 = true))"
            Base.depwarn(message, :init)
            Base.depwarn(message, :solve)
            aliases = ODEAliasSpecifier(alias_u0 = values(kwargs).alias_u0)
        else
            aliases = ODEAliasSpecifier(alias_u0 = nothing)
        end

        if haskey(kwargs, :alias_du0)
            message = "`alias_du0` keyword argument is deprecated, to set `alias_du0`,
            please use an ODEAliasSpecifier, e.g. `solve(prob, alias = ODEAliasSpecifier(alias_du0 = true))"
            Base.depwarn(message, :init)
            Base.depwarn(message, :solve)
            aliases = ODEAliasSpecifier(alias_u0 = aliases.alias_u0, alias_du0 = values(kwargs).alias_du0)
        else
            aliases = ODEAliasSpecifier(alias_u0 = aliases.alias_u0, alias_du0 = nothing)
        end
        
        aliases 

    else
         # If alias isa Bool, all fields of ODEAliases set to alias
        if alias isa Bool
            aliases = ODEAliasSpecifier(alias = alias)
        elseif alias isa ODEAliasSpecifier 
            aliases = alias
        end
    end

    if isnothing(aliases.alias_f) || aliases.alias_f 
        f = prob.f
    else
        f = deepcopy(prob.f)
    end

    if isnothing(aliases.alias_p) || aliases.alias_p 
        p = prob.p
    else
        p = recursivecopy(prob.p)
    end

    if !isnothing(aliases.alias_u0) && aliases.alias_u0
        u = prob.u0
    else
        u = recursivecopy(prob.u0)
    end

    if _alg isa DAEAlgorithm
        if !isnothing(aliases.alias_du0) && aliases.alias_du0
            du = prob.du0
        else
            du = recursivecopy(prob.du0)
        end
        duprev = recursivecopy(du)
    else
        du = nothing
        duprev = nothing
    end

    uType = typeof(u)
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    uEltypeNoUnits = recursive_unitless_eltype(u)
    tTypeNoUnits = typeof(one(tType))

    if prob isa DiscreteProblem
        abstol_internal = false
    elseif abstol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            abstol_internal = unitfulvalue(real(convert(uBottomEltype,
                oneunit(uBottomEltype) *
                1 // 10^6)))
        else
            abstol_internal = unitfulvalue.(real.(oneunit.(u) .* 1 // 10^6))
        end
    else
        abstol_internal = real.(abstol)
    end

    if prob isa DiscreteProblem
        reltol_internal = false
    elseif reltol === nothing
        if uBottomEltypeNoUnits == uBottomEltype
            reltol_internal = unitfulvalue(real(convert(uBottomEltype,
                oneunit(uBottomEltype) * 1 // 10^3)))
        else
            reltol_internal = unitfulvalue.(real.(oneunit.(u) .* 1 // 10^3))
        end
    else
        reltol_internal = real.(reltol)
    end

    dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
    # dtmin is all abs => does not care about sign already.

    if !isdae && isinplace(prob) && u isa AbstractArray && eltype(u) <: Number &&
       uBottomEltypeNoUnits == uBottomEltype && tType == tTypeNoUnits # Could this be more efficient for other arrays?
        rate_prototype = recursivecopy(u)
    elseif prob isa DAEProblem
        rate_prototype = prob.du0
    else
        if (uBottomEltypeNoUnits == uBottomEltype && tType == tTypeNoUnits) ||
           eltype(u) <: Enum
            rate_prototype = u
        else # has units!
            rate_prototype = u / oneunit(tType)
        end
    end
    rateType = typeof(rate_prototype) ## Can be different if united

    if isdae
        if uBottomEltype == uBottomEltypeNoUnits
            res_prototype = u
        else
            res_prototype = one(u)
        end
        resType = typeof(res_prototype)
    end

    if isnothing(aliases.alias_tstops) || aliases.alias_tstops
        tstops = tstops
    else
        tstops = recursivecopy(tstops)
    end

    if tstops isa AbstractArray || tstops isa Tuple || tstops isa Number
        _tstops = nothing
    else
        _tstops = tstops
        tstops = ()
    end
    tstops_internal = initialize_tstops(tType, tstops, d_discontinuities, tspan)
    saveat_internal = initialize_saveat(tType, saveat, tspan)
    d_discontinuities_internal = initialize_d_discontinuities(tType, d_discontinuities,
        tspan)

    callbacks_internal = CallbackSet(callback)

    max_len_cb = DiffEqBase.max_vector_callback_length_int(callbacks_internal)
    if max_len_cb !== nothing
        uBottomEltypeReal = real(uBottomEltype)
        if isinplace(prob)
            callback_cache = DiffEqBase.CallbackCache(u, max_len_cb, uBottomEltypeReal,
                uBottomEltypeReal)
        else
            callback_cache = DiffEqBase.CallbackCache(max_len_cb, uBottomEltypeReal,
                uBottomEltypeReal)
        end
    else
        callback_cache = nothing
    end

    ### Algorithm-specific defaults ###
    save_idxs, saved_subsystem = SciMLBase.get_save_idxs_and_saved_subsystem(prob, save_idxs)

    if save_idxs === nothing
        ksEltype = Vector{rateType}
    else
        ks_prototype = rate_prototype[save_idxs]
        ksEltype = Vector{typeof(ks_prototype)}
    end

    # Have to convert in case passed in wrong.
    if save_idxs === nothing
        timeseries = timeseries_init === () ? uType[] :
                     convert(Vector{uType}, timeseries_init)
    else
        u_initial = u[save_idxs]
        timeseries = timeseries_init === () ? typeof(u_initial)[] :
                     convert(Vector{uType}, timeseries_init)
    end

    ts = ts_init === () ? tType[] : convert(Vector{tType}, ts_init)
    ks = ks_init === () ? ksEltype[] : convert(Vector{ksEltype}, ks_init)
    alg_choice = _alg isa CompositeAlgorithm ? Int[] : nothing

    if (!adaptive || !isadaptive(_alg)) && save_everystep && tspan[2] - tspan[1] != Inf
        if dt == 0
            steps = length(tstops)
        else
            # For fixed dt, the only time dtmin makes sense is if it's smaller than eps().
            # Therefore user specified dtmin doesn't matter, but we need to ensure dt>=eps()
            # to prevent infinite loops.
            abs(dt) < dtmin &&
                throw(ArgumentError("Supplied dt is smaller than dtmin"))
            steps = ceil(Int, internalnorm((tspan[2] - tspan[1]) / dt, tspan[1]))
        end
        sizehint!(timeseries, steps + 1)
        sizehint!(ts, steps + 1)
        sizehint!(ks, steps + 1)
    elseif save_everystep
        sizehint!(timeseries, 50)
        sizehint!(ts, 50)
        sizehint!(ks, 50)
    elseif !isempty(saveat_internal)
        savelength = length(saveat_internal) + 1
        if save_start == false
            savelength -= 1
        end
        if save_end == false && prob.tspan[2] in saveat_internal.valtree
            savelength -= 1
        end
        sizehint!(timeseries, savelength)
        sizehint!(ts, savelength)
        sizehint!(ks, savelength)
    else
        sizehint!(timeseries, 2)
        sizehint!(ts, 2)
        sizehint!(ks, 2)
    end

    QT, EEstT = if tTypeNoUnits <: Integer
        typeof(qmin), typeof(qmin)
    elseif prob isa DiscreteProblem
        # The QT fields are not used for DiscreteProblems
        constvalue(tTypeNoUnits), constvalue(tTypeNoUnits)
    else
        typeof(DiffEqBase.value(internalnorm(u, t))), typeof(internalnorm(u, t))
    end

    k = rateType[]

    if uses_uprev(_alg, adaptive) || calck
        uprev = recursivecopy(u)
    else
        # Some algorithms do not use `uprev` explicitly. In that case, we can save
        # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
        uprev = u
    end
    if allow_extrapolation
        uprev2 = recursivecopy(u)
    else
        uprev2 = uprev
    end

    if prob isa DAEProblem
        cache = alg_cache(_alg, du, u, res_prototype, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt,
            reltol_internal, p, calck, Val(isinplace(prob)))
    else
        cache = alg_cache(_alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
            tTypeNoUnits, uprev, uprev2, f, t, dt, reltol_internal, p, calck,
            Val(isinplace(prob)))
    end

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
        controller = default_controller(_alg, cache, qoldinit, beta1, beta2)
    end

    save_end_user = save_end
    save_end = save_end === nothing ?
               save_everystep || isempty(saveat) || saveat isa Number ||
               prob.tspan[2] in saveat : save_end

    opts = DEOptions{typeof(abstol_internal), typeof(reltol_internal),
        QT, tType, typeof(controller),
        typeof(internalnorm), typeof(internalopnorm),
        typeof(save_end_user),
        typeof(callbacks_internal),
        typeof(isoutofdomain),
        typeof(progress_message), typeof(unstable_check),
        typeof(tstops_internal),
        typeof(d_discontinuities_internal), typeof(userdata),
        typeof(save_idxs),
        typeof(maxiters), typeof(tstops),
        typeof(saveat), typeof(d_discontinuities)}(maxiters, save_everystep,
        adaptive, abstol_internal,
        reltol_internal,
        QT(gamma), QT(qmax),
        QT(qmin),
        QT(qsteady_max),
        QT(qsteady_min),
        QT(qoldinit),
        QT(failfactor),
        tType(dtmax), tType(dtmin),
        controller,
        internalnorm,
        internalopnorm,
        save_idxs, tstops_internal,
        saveat_internal,
        d_discontinuities_internal,
        tstops, saveat,
        d_discontinuities,
        userdata, progress,
        progress_steps,
        progress_name,
        progress_message,
        progress_id,
        timeseries_errors,
        dense_errors, dense,
        save_on, save_start,
        save_end, save_end_user,
        callbacks_internal,
        isoutofdomain,
        unstable_check,
        verbose, calck, force_dtmin,
        advance_to_tstop,
        stop_at_next_tstop)

    stats = SciMLBase.DEStats(0)
    differential_vars = prob isa DAEProblem ? prob.differential_vars :
                        get_differential_vars(f, u)

    id = InterpolationData(
        f, timeseries, ts, ks, alg_choice, dense, cache, differential_vars, false)
    sol = DiffEqBase.build_solution(prob, _alg, ts, timeseries,
        dense = dense, k = ks, interp = id, alg_choice = alg_choice,
        calculate_error = false, stats = stats, saved_subsystem = saved_subsystem)

    if recompile_flag == true
        FType = typeof(f)
        SolType = typeof(sol)
        cacheType = typeof(cache)
    else
        FType = Function
        if _alg isa OrdinaryDiffEqAlgorithm
            SolType = DiffEqBase.AbstractODESolution
            cacheType = OrdinaryDiffEqCache
        else
            SolType = DiffEqBase.AbstractDAESolution
            cacheType = DAECache
        end
    end

    # rate/state = (state/time)/state = 1/t units, internalnorm drops units
    # we don't want to differentiate through eigenvalue estimation
    eigen_est = inv(one(tType))
    tprev = t
    dtcache = tType(dt)
    dtpropose = tType(dt)
    iter = 0
    kshortsize = 0
    reeval_fsal = false
    u_modified = false
    EEst = EEstT(1)
    just_hit_tstop = false
    isout = false
    accept_step = false
    force_stepfail = false
    last_stepfail = false
    do_error_check = true
    event_last_time = 0
    vector_event_last_time = 1
    last_event_error = prob isa DiscreteProblem ? false : zero(uBottomEltypeNoUnits)
    dtchangeable = isdtchangeable(_alg)
    q11 = QT(1)
    success_iter = 0
    erracc = QT(1)
    dtacc = tType(1)
    reinitiailize = true
    saveiter = 0 # Starts at 0 so first save is at 1
    saveiter_dense = 0
    fsalfirst, fsallast = get_fsalfirstlast(cache, rate_prototype)

    integrator = ODEIntegrator{typeof(_alg), isinplace(prob), uType, typeof(du),
        tType, typeof(p),
        typeof(eigen_est), typeof(EEst),
        QT, typeof(tdir), typeof(k), SolType,
        FType, cacheType,
        typeof(opts), typeof(fsalfirst),
        typeof(last_event_error), typeof(callback_cache),
        typeof(initializealg), typeof(differential_vars)}(
        sol, u, du, k, t, tType(dt), f, p,
        uprev, uprev2, duprev, tprev,
        _alg, dtcache, dtchangeable,
        dtpropose, tdir, eigen_est, EEst,
        QT(qoldinit), q11,
        erracc, dtacc, success_iter,
        iter, saveiter, saveiter_dense, cache,
        callback_cache,
        kshortsize, force_stepfail,
        last_stepfail,
        just_hit_tstop, do_error_check,
        event_last_time,
        vector_event_last_time,
        last_event_error, accept_step,
        isout, reeval_fsal,
        u_modified, reinitiailize, isdae,
        opts, stats, initializealg, differential_vars,
        fsalfirst, fsallast)

    if initialize_integrator
        if isdae || SciMLBase.has_initializeprob(prob.f)
            DiffEqBase.initialize_dae!(integrator)
            update_uprev!(integrator)
        end

        if save_start
            integrator.saveiter += 1 # Starts at 1 so first save is at 2
            integrator.saveiter_dense += 1
            copyat_or_push!(ts, 1, t)
            # N.B.: integrator.u can be modified by initialized_dae!
            if save_idxs === nothing
                copyat_or_push!(timeseries, 1, integrator.u)
                copyat_or_push!(ks, 1, [rate_prototype])
            else
                copyat_or_push!(timeseries, 1, integrator.u[save_idxs], Val{false})
                copyat_or_push!(ks, 1, [ks_prototype])
            end
        else
            integrator.saveiter = 0 # Starts at 0 so first save is at 1
            integrator.saveiter_dense = 0
        end

        initialize_callbacks!(integrator, initialize_save)
        initialize!(integrator, integrator.cache)

        if _alg isa OrdinaryDiffEqCompositeAlgorithm
            # in case user mixes adaptive and non-adaptive algorithms
            ensure_behaving_adaptivity!(integrator, integrator.cache)

            if save_start
                # Loop to get all of the extra possible saves in callback initialization
                for i in 1:(integrator.saveiter)
                    copyat_or_push!(alg_choice, i, integrator.cache.current)
                end
            end
        end
    end

    if _tstops !== nothing
        tstops = _tstops(parameter_values(integrator), prob.tspan)
        for tstop in tstops
            add_tstop!(integrator, tstop)
        end
    end

    handle_dt!(integrator)
    integrator
end

function DiffEqBase.solve!(integrator::ODEIntegrator)
    @inbounds while !isempty(integrator.opts.tstops)
        while integrator.tdir * integrator.t < first(integrator.opts.tstops)
            loopheader!(integrator)
            if integrator.do_error_check && check_error!(integrator) != ReturnCode.Success
                return integrator.sol
            end
            perform_step!(integrator, integrator.cache)
            loopfooter!(integrator)
            if isempty(integrator.opts.tstops)
                break
            end
        end
        handle_tstop!(integrator)
    end
    postamble!(integrator)

    f = integrator.sol.prob.f

    if DiffEqBase.has_analytic(f)
        DiffEqBase.calculate_solution_errors!(integrator.sol;
            timeseries_errors = integrator.opts.timeseries_errors,
            dense_errors = integrator.opts.dense_errors)
    end
    if integrator.sol.retcode != ReturnCode.Default
        return integrator.sol
    end
    integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, ReturnCode.Success)
end

# Helpers

function handle_dt!(integrator)
    if iszero(integrator.dt) && integrator.opts.adaptive
        auto_dt_reset!(integrator)
        if sign(integrator.dt) != integrator.tdir && !iszero(integrator.dt) &&
           !isnan(integrator.dt)
            error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
        end
        if isnan(integrator.dt)
            if integrator.opts.verbose
                @warn("Automatic dt set the starting dt as NaN, causing instability. Exiting.")
            end
        end
    elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) &&
           integrator.tdir < 0
        integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
    end
end

# time stops
@inline function initialize_tstops(::Type{T}, tstops, d_discontinuities, tspan) where {T}
    tstops_internal = BinaryHeap{T}(DataStructures.FasterForward())

    t0, tf = tspan
    tdir = sign(tf - t0)
    tdir_t0 = tdir * t0
    tdir_tf = tdir * tf

    for t in tstops
        tdir_t = tdir * t
        tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end
    for t in d_discontinuities
        tdir_t = tdir * t
        tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end
    push!(tstops_internal, tdir_tf)

    return tstops_internal
end

# saving time points
function initialize_saveat(::Type{T}, saveat, tspan) where {T}
    saveat_internal = BinaryHeap{T}(DataStructures.FasterForward())

    t0, tf = tspan
    tdir = sign(tf - t0)
    tdir_t0 = tdir * t0
    tdir_tf = tdir * tf

    if saveat isa Number
        directional_saveat = tdir * abs(saveat)
        for t in (t0 + directional_saveat):directional_saveat:tf
            push!(saveat_internal, tdir * t)
        end
    elseif !isempty(saveat)
        for t in saveat
            tdir_t = tdir * t
            tdir_t0 < tdir_t ≤ tdir_tf && push!(saveat_internal, tdir_t)
        end
    end

    return saveat_internal
end

# discontinuities
function initialize_d_discontinuities(::Type{T}, d_discontinuities, tspan) where {T}
    d_discontinuities_internal = BinaryHeap{T}(DataStructures.FasterForward())
    sizehint!(d_discontinuities_internal, length(d_discontinuities))

    t0, tf = tspan
    tdir = sign(tf - t0)

    for t in d_discontinuities
        push!(d_discontinuities_internal, tdir * t)
    end

    return d_discontinuities_internal
end

function initialize_callbacks!(integrator, initialize_save = true)
    t = integrator.t
    u = integrator.u
    callbacks = integrator.opts.callback
    integrator.u_modified = true

    u_modified = initialize!(callbacks, u, t, integrator)

    # if the user modifies u, we need to fix previous values before initializing
    # FSAL in order for the starting derivatives to be correct
    if u_modified
        if isinplace(integrator.sol.prob)
            recursivecopy!(integrator.uprev, integrator.u)
        else
            integrator.uprev = integrator.u
        end

        if alg_extrapolates(integrator.alg)
            if isinplace(integrator.sol.prob)
                recursivecopy!(integrator.uprev2, integrator.uprev)
            else
                integrator.uprev2 = integrator.uprev
            end
        end

        if initialize_save &&
           (any((c) -> c.save_positions[2], callbacks.discrete_callbacks) ||
            any((c) -> c.save_positions[2], callbacks.continuous_callbacks))
            savevalues!(integrator, true)
        end
    end

    # reset this as it is now handled so the integrators should proceed as normal
    integrator.u_modified = false
end

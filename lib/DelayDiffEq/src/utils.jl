"""
    has_constant_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains constant delays.
"""
has_constant_lags(integrator::DDEIntegrator) = has_constant_lags(integrator.sol.prob)

"""
    has_dependent_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains dependent delays.
"""
has_dependent_lags(integrator::DDEIntegrator) = has_dependent_lags(integrator.sol.prob)

"""
    has_constant_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains constant delays.
"""
function has_constant_lags(prob::DDEProblem)
    return prob.constant_lags !== nothing && !isempty(prob.constant_lags)
end

"""
    has_dependent_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains dependent delays.
"""
function has_dependent_lags(prob::DDEProblem)
    return prob.dependent_lags !== nothing && !isempty(prob.dependent_lags)
end

"""
    u_uprev(u0, alg; kwargs...)

Return state vectors `u` and `uprev` (possibly aliased) for solving the
differential equation problem for initial state `u0` with algorithm `alg`.
"""
function u_uprev(
        u0, alg;
        alias_u0 = false,
        adaptive = isadaptive(alg),
        calck = false
    )
    if alias_u0
        u = u0
    else
        u = recursivecopy(u0)
    end

    # Some algorithms do not use `uprev` explicitly. In that case, we can save
    # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
    if OrdinaryDiffEqCore.uses_uprev(alg, adaptive) || calck
        uprev = recursivecopy(u)
    else
        uprev = u
    end

    return u, uprev
end

"""
    u_uprev_uprev2(u0, alg; kwargs...)

Return state vectors `u`, `uprev`, and `uprev2` (possibly aliased) for solving the
differential equation problem for initial state `u0` with algorithm `alg`.
"""
function u_uprev_uprev2(
        u0, alg;
        allow_extrapolation = alg_extrapolates(alg),
        kwargs...
    )
    # compute u and uprev first
    u, uprev = u_uprev(u0, alg; kwargs...)

    if allow_extrapolation
        uprev2 = recursivecopy(u)
    else
        uprev2 = uprev
    end

    return u, uprev, uprev2
end

"""
    get_abstol(u, tspan, alg; abstol = nothing)

Return the absolute tolerance for solving the differential equation problem with state
variable `u` and time span `tspan` with algorithm `alg`.
"""
function get_abstol(u, tspan, alg; abstol = nothing)
    if alg isa FunctionMap
        _abstol = real.(zero.(u))
    elseif abstol === nothing
        uBottomEltype = recursive_bottom_eltype(u)
        uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

        if uBottomEltypeNoUnits == uBottomEltype
            _abstol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^6))
        else
            _abstol = real.(oneunit.(u) .* 1 // 10^6)
        end
    else
        _abstol = real.(abstol)
    end

    return _abstol
end

"""
    get_reltol(u, tspan, alg; reltol = nothing)

Return the relative tolerance for solving the differential equation problem with state
variable `u` and time span `tspan` with algorithm `alg`.
"""
function get_reltol(u, tspan, alg; reltol = nothing)
    if alg isa FunctionMap
        _reltol = real.(zero(first(u) / t))
    elseif reltol === nothing
        uBottomEltype = recursive_bottom_eltype(u)
        uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

        if uBottomEltypeNoUnits == uBottomEltype
            _reltol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^3))
        else
            _reltol = real.(oneunit.(u) .* 1 // 10^3)
        end
    else
        _reltol = real.(reltol)
    end

    return _reltol
end

"""
    callback_set_and_cache(prob, callback)

Return set of callbacks and its cache for the differential equation problem `prob` and the
user-provided `callback`.
"""
function callback_set_and_cache(prob, callback)
    callback_set = CallbackSet(callback)

    max_len_cb = DiffEqBase.max_vector_callback_length(callback_set)
    if max_len_cb isa VectorContinuousCallback
        uBottomEltype = recursive_bottom_eltype(prob.u0)
        callback_cache = DiffEqBase.CallbackCache(
            max_len_cb.len, uBottomEltype,
            uBottomEltype
        )
    else
        callback_cache = nothing
    end

    return callback_set, callback_cache
end

"""
    rate_prototype_of(u0, tspan)

Return prototype of rates for a given differential equation problem with state `u` and
time span `tspan`.
"""
rate_prototype_of(u0, tspan) = @.. u0 * $(inv(oneunit(eltype(tspan))))

"""
    solution_arrays(u, tspan, rate_prototype; kwargs...)

Return arrays of saved time points, states, and rates, initialized with the solution at the
first time point if `save_start = true` (the default).
"""
function solution_arrays(
        u, tspan, rate_prototype;
        timeseries_init,
        ts_init,
        ks_init,
        save_idxs,
        save_start
    )
    # determine types of time and state
    uType = typeof(u)
    tType = eltype(tspan)

    # initialize vector of saved time points
    ts = ts_init === () ? tType[] : convert(Vector{tType}, ts_init)

    # initialize vector of saved states
    if save_idxs === nothing
        timeseries = timeseries_init === () ? uType[] :
            convert(Vector{uType}, timeseries_init)
    else
        u_initial = u[save_idxs]
        timeseries = timeseries_init === () ? typeof(u_initial)[] :
            convert(Vector{typeof(u_initial)}, timeseries_init)
    end

    # initialize vector of saved rates
    if save_idxs === nothing
        ksEltype = Vector{typeof(rate_prototype)}
    else
        ks_prototype = rate_prototype[save_idxs]
        ksEltype = Vector{typeof(ks_prototype)}
    end
    ks = ks_init === () ? ksEltype[] : convert(Vector{ksEltype}, ks_init)

    # save solution at initial time point
    if save_start
        copyat_or_push!(ts, 1, first(tspan))
        if save_idxs === nothing
            copyat_or_push!(timeseries, 1, u)
            copyat_or_push!(ks, 1, [rate_prototype])
        else
            u_initial = u[save_idxs]
            copyat_or_push!(timeseries, 1, u_initial, Val{false})
            copyat_or_push!(ks, 1, [ks_prototype])
        end
    end

    return ts, timeseries, ks
end

"""
    _sizehint_solution!(sol::DESolution, n)

Suggest that solution `sol` reserves capacity for at least `n` elements.
"""
function _sizehint_solution!(sol::DESolution, n)
    sizehint!(sol.u, n)
    sizehint!(sol.t, n)
    sizehint!(sol.k, n)

    return nothing
end

"""
    _sizehint_solution!(sol::DESolution, alg, tspan, tstops, saveat; kwargs...)

Suggest that solution `sol` reserves capacity for a number of elements that
depends on the parameter settings of the numerical solver.
"""
function _sizehint_solution!(
        sol::DESolution, alg, tspan, tstops, saveat;
        save_everystep, adaptive, internalnorm, dt, dtmin
    )
    # obtain integration time
    t0 = first(tspan)
    integrationtime = last(tspan) - t0

    if !adaptive && save_everystep && !isinf(integrationtime)
        # determine number of steps if known a priori
        if iszero(dt)
            steps = length(tstops)
        else
            abs(dt) < dtmin && throw(ArgumentError("Supplied dt is smaller than dtmin"))
            steps = ceil(Int, internalnorm(integrationtime / dt, t0))
        end
        _sizehint_solution!(sol, steps + 1)
    elseif save_everystep
        _sizehint_solution!(sol, 50)
    elseif !isempty(saveat)
        _sizehint_solution!(sol, length(saveat) + 1)
    else
        _sizehint_solution!(sol, 2)
    end

    return nothing
end

function build_history_function(
        prob, alg, rate_prototype, reltol, differential_vars;
        dt, dtmin, adaptive, calck, internalnorm
    )
    (; f, u0, tspan, p) = prob

    t0 = first(tspan)
    tType = eltype(tspan)
    tTypeNoUnits = typeof(one(tType))
    tdir = sign(last(tspan) - t0)

    uEltypeNoUnits = recursive_unitless_eltype(u0)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u0)

    # bootstrap an ODE integrator
    # - whose solution captures the dense history of the simulation
    # - that is used for extrapolation of the history for time points past the
    #   already fixed history
    # - that is used for interpolation of the history for time points in the
    #   current integration step (so the interpolation is fixed while updating the stages)
    # we wrap the user-provided history function such that function calls during the setup
    # of the integrator do not fail
    ode_f = ODEFunctionWrapper(f, prob.h)
    ode_prob = ODEProblem{isinplace(prob)}(ode_f, u0, tspan, p)

    # get states of ODE integrator (do not alias uprev)
    ode_u, ode_uprev = u_uprev(u0, alg; alias_u0 = false, calck = true)

    # initialize output arrays
    ode_k = typeof(rate_prototype)[]
    ode_ts, ode_timeseries,
        ode_ks = solution_arrays(
        ode_u, tspan, rate_prototype;
        timeseries_init = (),
        ts_init = (),
        ks_init = (),
        save_idxs = nothing,
        save_start = true
    )

    # obtain cache (we alias uprev2 and uprev)
    ode_cache = OrdinaryDiffEqCore.alg_cache(
        alg.alg, ode_u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, ode_uprev,
        ode_uprev, ode_f, t0, zero(tType), reltol, p,
        calck,
        Val(isinplace(prob)), OrdinaryDiffEqCore.DEVerbosity()
    )

    # build dense interpolation of history
    ode_alg_choice = iscomposite(alg) ? Int[] : nothing
    ode_id = OrdinaryDiffEqCore.InterpolationData(
        ode_f, ode_timeseries, ode_ts,
        ode_ks,
        ode_alg_choice, true, ode_cache,
        differential_vars, false
    )
    ode_sol = SciMLBase.build_solution(
        ode_prob, alg.alg, ode_ts, ode_timeseries;
        dense = true, k = ode_ks, interp = ode_id,
        alg_choice = ode_alg_choice,
        calculate_error = false,
        stats = DiffEqBase.Stats(0)
    )

    # reserve capacity
    _sizehint_solution!(
        ode_sol, alg.alg, tspan, (), ();
        save_everystep = true, adaptive = adaptive, internalnorm = internalnorm,
        dt = dt, dtmin = dtmin
    )

    # create simple integrator
    tdirType = typeof(sign(zero(tType)))
    ode_integrator = HistoryODEIntegrator{
        typeof(alg.alg), isinplace(prob), typeof(prob.u0),
        tType, tdirType, typeof(ode_k),
        typeof(ode_sol), typeof(ode_cache),
        typeof(differential_vars),
    }(
        ode_sol,
        ode_u, ode_k,
        t0,
        zero(tType),
        ode_uprev,
        t0, alg.alg,
        zero(tType),
        tdir, 1, 1,
        ode_cache,
        differential_vars
    )

    # combine the user-provided history function and the ODE integrator with dense solution
    # to a joint dense history of the DDE
    # we use this history information to create a problem function of the DDE with all
    # available history information that is of the form f(du,u,p,t) or f(u,p,t) such that
    # ODE algorithms can be applied
    return HistoryFunction(prob.h, ode_integrator)
end

"""
    initialize_solution!(integrator::DDEIntegrator)

Initialize the solution of an integrator by adjusting the cache for composite algorithms.
"""
function initialize_solution!(integrator::DDEIntegrator)
    if iscomposite(integrator.alg)
        copyat_or_push!(integrator.integrator.sol.alg_choice, 1, integrator.cache.current)
        if integrator.opts.save_start
            copyat_or_push!(integrator.sol.alg_choice, 1, integrator.cache.current)
        end
    end

    return nothing
end

function unwrap_alg(integrator::DDEIntegrator, is_stiff)
    alg = integrator.alg
    iscomp = alg isa CompositeAlgorithm
    if !iscomp
        return alg
    elseif alg.choice_function isa AutoSwitch
        num = is_stiff ? 2 : 1
        return alg.algs[num]
    else
        return alg.algs[integrator.cache.current]
    end
end

function OrdinaryDiffEqCore.nlsolve_f(integrator::DDEIntegrator)
    return OrdinaryDiffEqCore.nlsolve_f(integrator.f, unwrap_alg(integrator, true))
end

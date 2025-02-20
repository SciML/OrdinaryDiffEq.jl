# We want to make sure that the first argument of change_t_via_interpolation!
# is specialized, yet, it needs to take both ODEIntegrator and DDEIntegrator.
# Hence, we need to have two separate functions.

function _change_t_via_interpolation!(integrator, t,
        modify_save_endpoint::Type{Val{T}}, reinitialize_alg = nothing) where {T}
    # Can get rid of an allocation here with a function
    # get_tmp_arr(integrator.cache) which gives a pointer to some
    # cache array which can be modified.
    if integrator.tdir * t < integrator.tdir * integrator.tprev
        error("Current interpolant only works between tprev and t")
    elseif t != integrator.t
        if is_constant_cache(integrator.cache)
            integrator.u = integrator(t)
        else
            integrator(integrator.u, t)
        end
        integrator.t = t
        integrator.dt = integrator.t - integrator.tprev
        DiffEqBase.reeval_internals_due_to_modification!(
            integrator; callback_initializealg = reinitialize_alg)
        if T
            solution_endpoint_match_cur_integrator!(integrator)
        end
    end
    return nothing
end
function DiffEqBase.change_t_via_interpolation!(integrator::ODEIntegrator,
        t,
        modify_save_endpoint::Type{Val{T}} = Val{
            false,
        }, reinitialize_alg = nothing) where {
        T,
}
    _change_t_via_interpolation!(integrator, t, modify_save_endpoint, reinitialize_alg)
    return nothing
end

function DiffEqBase.reeval_internals_due_to_modification!(
        integrator::ODEIntegrator, continuous_modification = true;
        callback_initializealg = nothing)
    if integrator.isdae
        DiffEqBase.initialize_dae!(integrator,
            isnothing(callback_initializealg) ? integrator.initializealg :
            callback_initializealg)
        update_uprev!(integrator)
    end

    if continuous_modification && integrator.opts.calck
        resize!(integrator.k, integrator.kshortsize) # Reset k for next step!
        alg = unwrap_alg(integrator, false)
        if SciMLBase.has_lazy_interpolation(alg)
            ode_addsteps!(integrator, integrator.f, true, false, !alg.lazy)
        else
            ode_addsteps!(integrator, integrator.f, true, false)
        end
    end

    integrator.u_modified = false
    integrator.reeval_fsal = true
end

@inline function DiffEqBase.get_du(integrator::ODEIntegrator)
    isdiscretecache(integrator.cache) &&
        error("Derivatives are not defined for this stepper.")
    return if isfsal(integrator.alg)
        integrator.fsallast
    else
        integrator(integrator.t, Val{1})
    end
end

@inline function DiffEqBase.get_du!(out, integrator::ODEIntegrator)
    isdiscretecache(integrator.cache) &&
        error("Derivatives are not defined for this stepper.")
    if isdiscretecache(integrator.cache)
        out .= integrator.cache.tmp
    else
        return if isfsal(integrator.alg) &&
                  !has_stiff_interpolation(integrator.alg)
            # Special stiff interpolations do not store the
            # right value in fsallast
            out .= integrator.fsallast
        else
            integrator(out, integrator.t, Val{1})
        end
    end
end

function u_modified!(integrator::ODEIntegrator, bool::Bool)
    integrator.u_modified = bool
end

function get_proposed_dt(integrator::ODEIntegrator)
    ifelse(integrator.opts.adaptive, integrator.dtpropose, integrator.dtcache)
end
function set_proposed_dt!(integrator::ODEIntegrator, dt::Number)
    (integrator.dtpropose = dt; integrator.dtcache = dt)
end

function set_proposed_dt!(integrator::ODEIntegrator, integrator2::ODEIntegrator)
    integrator.dtpropose = integrator2.dtpropose
    integrator.dtcache = integrator2.dtcache
    integrator.qold = integrator2.qold
    integrator.erracc = integrator2.erracc
    integrator.dtacc = integrator2.dtacc
end

#TODO: Bigger caches for most algorithms
@inline function DiffEqBase.get_tmp_cache(integrator::ODEIntegrator)
    get_tmp_cache(integrator::ODEIntegrator, integrator.alg, integrator.cache)
end

# the ordering of the cache arrays is important!!!
@inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqAlgorithm,
        cache::OrdinaryDiffEqConstantCache)
    nothing
end
@inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp,)
end
@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.nlsolver.tmp, cache.atmp)
end
@inline function DiffEqBase.get_tmp_cache(integrator, alg::OrdinaryDiffEqNewtonAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.nlsolver.tmp, cache.nlsolver.z)
end
@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp, cache.linsolve_tmp)
end

@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp, cache.utilde)
end
@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp, cache.dz)
end
@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqLinearExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp,)
end
@inline function DiffEqBase.get_tmp_cache(integrator, alg::CompositeAlgorithm,
        cache::CompositeCache)
    get_tmp_cache(integrator, alg.algs[1], cache.caches[1])
end
@inline function DiffEqBase.get_tmp_cache(integrator, alg::CompositeAlgorithm,
        cache::DefaultCache)
    init_ith_default_cache(cache, alg.algs, cache.current)
    if cache.current == 1
        get_tmp_cache(integrator, alg.algs[1], cache.cache1)
    elseif cache.current == 2
        get_tmp_cache(integrator, alg.algs[2], cache.cache2)
    elseif cache.current == 3
        get_tmp_cache(integrator, alg.algs[3], cache.cache3)
    elseif cache.current == 4
        get_tmp_cache(integrator, alg.algs[4], cache.cache4)
    elseif cache.current == 5
        get_tmp_cache(integrator, alg.algs[5], cache.cache5)
    else
        @assert cache.current == 6
        get_tmp_cache(integrator, alg.algs[6], cache.cache6)
    end
end

@inline function DiffEqBase.get_tmp_cache(integrator, alg::DAEAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.nlsolver.cache.dz, cache.atmp)
end

function full_cache(integrator::ODEIntegrator)
    # for DefaultCache, we need to make sure to initialize all the caches in case they get switched to later
    if integrator.cache isa DefaultCache
        @unpack alg, cache = integrator
        algs = alg.algs
        init_ith_default_cache(cache, algs, 1)
        init_ith_default_cache(cache, algs, 2)
        init_ith_default_cache(cache, algs, 3)
        init_ith_default_cache(cache, algs, 4)
        init_ith_default_cache(cache, algs, 5)
        init_ith_default_cache(cache, algs, 6)
    end
    full_cache(integrator.cache)
end
function full_cache(cache::CompositeCache)
    Iterators.flatten(full_cache(c) for c in cache.caches)
end
function full_cache(cache::DefaultCache)
    caches = (
        cache.cache1, cache.cache2, cache.cache3, cache.cache4, cache.cache5, cache.cache6)
    Iterators.flatten(full_cache(c) for c in caches)
end

function DiffEqBase.add_tstop!(integrator::ODEIntegrator, t)
    integrator.tdir * (t - integrator.t) < zero(integrator.t) &&
        error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
    push!(integrator.opts.tstops, integrator.tdir * t)
end

DiffEqBase.has_tstop(integrator::ODEIntegrator) = !isempty(integrator.opts.tstops)
DiffEqBase.first_tstop(integrator::ODEIntegrator) = first(integrator.opts.tstops)
DiffEqBase.pop_tstop!(integrator::ODEIntegrator) = pop!(integrator.opts.tstops)

function DiffEqBase.add_saveat!(integrator::ODEIntegrator, t)
    integrator.tdir * (t - integrator.t) < zero(integrator.t) &&
        error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
    push!(integrator.opts.saveat, integrator.tdir * t)
end

function resize!(integrator::ODEIntegrator, i::Int)
    @unpack cache = integrator

    for c in full_cache(integrator)
        # Skip nothings which may exist in the cache since extra variables
        # may be required for things like units
        c !== nothing && resize!(c, i)
    end
    !isnothing(integrator.fsalfirst) && resize!(integrator.fsalfirst, i)
    !isnothing(integrator.fsallast) && resize!(integrator.fsallast, i)
    resize_f!(integrator.f, i)
    resize_nlsolver!(integrator, i)
    resize_J_W!(cache, integrator, i)
    resize_non_user_cache!(integrator, cache, i)
end
# we can't use resize!(..., i::Union{Int, NTuple{N,Int}}) where {N} because of method ambiguities with DiffEqBase
function resize!(integrator::ODEIntegrator, i::NTuple{N, Int}) where {N}
    @unpack cache = integrator

    for c in full_cache(cache)
        resize!(c, i)
    end
    !isnothing(integrator.fsalfirst) && resize!(integrator.fsalfirst, i)
    !isnothing(integrator.fsallast) && resize!(integrator.fsallast, i)
    resize_f!(integrator.f, i)
    # TODO the parts below need to be adapted for implicit methods
    isdefined(integrator.cache, :nlsolver) && resize_nlsolver!(integrator, i)
    resize_J_W!(cache, integrator, i)
    resize_non_user_cache!(integrator, cache, i)
end

# default fallback
resize_f!(f, i) = nothing

function resize_f!(f::SplitFunction, i)
    resize!(f._func_cache, i)
    return nothing
end

function resize_J_W! end

function resize_non_user_cache!(integrator::ODEIntegrator, i::Int)
    resize_non_user_cache!(integrator, integrator.cache, i)
end
function deleteat_non_user_cache!(integrator::ODEIntegrator, i)
    deleteat_non_user_cache!(integrator, integrator.cache, i)
end
function addat_non_user_cache!(integrator::ODEIntegrator, i)
    addat_non_user_cache!(integrator, integrator.cache, i)
end

resize_non_user_cache!(integrator::ODEIntegrator, cache, i) = nothing

function resize_non_user_cache!(integrator::ODEIntegrator, cache::CompositeCache, i)
    for _cache in cache.caches
        resize_non_user_cache!(integrator, _cache, i)
    end
end

function deleteat_non_user_cache!(integrator::ODEIntegrator, cache::CompositeCache, i)
    for _cache in cache.caches
        deleteat_non_user_cache!(integrator, _cache, i)
    end
end

function addat_non_user_cache!(integrator::ODEIntegrator, cache::CompositeCache, i)
    for _cache in cache.caches
        addat_non_user_cache!(integrator, _cache, i)
    end
end

function deleteat_non_user_cache!(integrator::ODEIntegrator, cache, idxs)
    # ordering doesn't matter in deterministic cache, so just resize
    # to match the size of u
    i = length(integrator.u)
    resize_non_user_cache!(integrator, cache, i)
end

function addat_non_user_cache!(integrator::ODEIntegrator, cache, idxs)
    # ordering doesn't matter in deterministic cache, so just resize
    # to match the size of u
    i = length(integrator.u)
    resize_non_user_cache!(integrator, cache, i)
end

function deleteat!(integrator::ODEIntegrator, idxs)
    for c in full_cache(integrator)
        deleteat!(c, idxs)
    end
    deleteat_non_user_cache!(integrator, integrator.cache, idxs)
end

function addat!(integrator::ODEIntegrator, idxs)
    for c in full_cache(integrator)
        addat!(c, idxs)
    end
    addat_non_user_cache!(integrator, integrator.cache, idxs)
end

function terminate!(integrator::ODEIntegrator, retcode = ReturnCode.Terminated)
    integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, retcode)
    integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end

DiffEqBase.has_reinit(integrator::ODEIntegrator) = true
function DiffEqBase.reinit!(integrator::ODEIntegrator, u0 = integrator.sol.prob.u0;
        t0 = integrator.sol.prob.tspan[1],
        tf = integrator.sol.prob.tspan[2],
        erase_sol = true,
        tstops = integrator.opts.tstops_cache,
        saveat = integrator.opts.saveat_cache,
        d_discontinuities = integrator.opts.d_discontinuities_cache,
        reset_dt = (integrator.dtcache == zero(integrator.dt)) &&
            integrator.opts.adaptive,
        reinit_dae = true,
        reinit_callbacks = true, initialize_save = true,
        reinit_cache = true,
        reinit_retcode = true)
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.u, u0)
        recursivecopy!(integrator.uprev, integrator.u)
    else
        integrator.u = u0
        integrator.uprev = integrator.u
    end

    if alg_extrapolates(integrator.alg)
        if isinplace(integrator.sol.prob)
            recursivecopy!(integrator.uprev2, integrator.uprev)
        else
            integrator.uprev2 = integrator.uprev
        end
    end

    integrator.t = t0
    integrator.tprev = t0

    tType = typeof(integrator.t)
    tspan = (tType(t0), tType(tf))
    integrator.opts.tstops = initialize_tstops(tType, tstops, d_discontinuities, tspan)
    integrator.opts.saveat = initialize_saveat(tType, saveat, tspan)
    integrator.opts.d_discontinuities = initialize_d_discontinuities(tType,
        d_discontinuities,
        tspan)

    if erase_sol
        if integrator.opts.save_start
            resize_start = 1
        else
            resize_start = 0
        end
        resize!(integrator.sol.u, resize_start)
        resize!(integrator.sol.t, resize_start)
        resize!(integrator.sol.k, resize_start)

        if integrator.opts.save_start || (!isempty(saveat) && saveat[1] == tType(t0))
            copyat_or_push!(integrator.sol.t, 1, t0)
            if integrator.opts.save_idxs === nothing
                copyat_or_push!(integrator.sol.u, 1, u0)
            else
                u_initial = u0[integrator.opts.save_idxs]
                copyat_or_push!(integrator.sol.u, 1, u_initial, Val{false})
            end
        end
        if integrator.sol.u_analytic !== nothing
            resize!(integrator.sol.u_analytic, 0)
        end
        if integrator.alg isa OrdinaryDiffEqCompositeAlgorithm
            resize!(integrator.sol.alg_choice, resize_start)
        end
        integrator.saveiter = resize_start
        if integrator.opts.dense
            integrator.saveiter_dense = resize_start
        end
    end
    integrator.iter = 0
    integrator.success_iter = 0
    integrator.u_modified = false

    # full re-initialize the PI in timestepping
    reinit!(integrator, integrator.opts.controller)
    integrator.qold = integrator.opts.qoldinit
    integrator.q11 = typeof(integrator.q11)(1)
    integrator.erracc = typeof(integrator.erracc)(1)
    integrator.dtacc = typeof(integrator.dtacc)(1)

    if reset_dt
        auto_dt_reset!(integrator)
    end

    if reinit_dae &&
       (integrator.isdae || SciMLBase.has_initializeprob(integrator.sol.prob.f))
        DiffEqBase.initialize_dae!(integrator)
        update_uprev!(integrator)
    end

    if reinit_callbacks
        initialize_callbacks!(integrator, initialize_save)
    end

    if reinit_cache
        initialize!(integrator, integrator.cache)
    end

    if reinit_retcode
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, ReturnCode.Default)
    end
    return nothing
end

function DiffEqBase.auto_dt_reset!(integrator::ODEIntegrator)
    integrator.dt = ode_determine_initdt(integrator.u, integrator.t,
        integrator.tdir, integrator.opts.dtmax,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, integrator.sol.prob,
        integrator)
    integrator.dtpropose = integrator.dt
    increment_nf!(integrator.stats, 2)
end

function increment_nf!(stats, amt = 1)
    stats.nf += amt
end

function DiffEqBase.set_t!(integrator::ODEIntegrator, t::Real)
    if integrator.opts.save_everystep
        error("Integrator time cannot be reset unless it is initialized",
            " with save_everystep=false")
    end
    if alg_extrapolates(integrator.alg) || !isdtchangeable(integrator.alg)
        reinit!(integrator, integrator.u;
            t0 = t,
            reset_dt = false,
            reinit_callbacks = false,
            reinit_cache = false)
    else
        integrator.t = t
    end
end

function DiffEqBase.set_u!(integrator::ODEIntegrator, u)
    if integrator.opts.save_everystep
        error("Integrator state cannot be reset unless it is initialized",
            " with save_everystep=false")
    end
    integrator.u = u
    u_modified!(integrator, true)
end

DiffEqBase.has_stats(i::ODEIntegrator) = true

DiffEqBase.get_tstops(integ::ODEIntegrator) = integ.opts.tstops
DiffEqBase.get_tstops_array(integ::ODEIntegrator) = get_tstops(integ).valtree
DiffEqBase.get_tstops_max(integ::ODEIntegrator) = maximum(get_tstops_array(integ))

# We want to make sure that the first argument of change_t_via_interpolation!
# is specialized, yet, it needs to take both ODEIntegrator and DDEIntegrator.
# Hence, we need to have two separate functions.

function _change_t_via_interpolation!(
        integrator, t,
        modify_save_endpoint::Type{Val{T}}, reinitialize_alg = nothing
    ) where {T}
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
        SciMLBase.reeval_internals_due_to_modification!(
            integrator; callback_initializealg = reinitialize_alg
        )
        if T
            solution_endpoint_match_cur_integrator!(integrator)
        end
    end
    return nothing
end
function SciMLBase.change_t_via_interpolation!(
        integrator::ODEIntegratorType,
        t,
        modify_save_endpoint::Type{Val{T}} = Val{
            false,
        }, reinitialize_alg = nothing
    ) where {
        T,
    }
    _change_t_via_interpolation!(integrator, t, modify_save_endpoint, reinitialize_alg)
    return nothing
end

function SciMLBase.reeval_internals_due_to_modification!(
        integrator::ODEIntegratorType, continuous_modification = true;
        callback_initializealg = nothing
    )
    if integrator.isdae
        DiffEqBase.initialize_dae!(
            integrator,
            isnothing(callback_initializealg) ? integrator.initializealg :
                callback_initializealg
        )
        update_uprev!(integrator)
    end

    if continuous_modification && integrator.opts.calck
        resize!(integrator.k, integrator.kshortsize) # Reset k for next step!
        alg = unwrap_alg(integrator, false)
        if SciMLBase.has_lazy_interpolation(alg)
            ode_addsteps!(integrator, integrator.f, true, false, !_unwrap_val(alg.lazy))
        else
            ode_addsteps!(integrator, integrator.f, true, false)
        end
    end

    integrator.u_modified = false
    return integrator.reeval_fsal = true
end

@inline function SciMLBase.get_du(integrator::ODEIntegratorType)
    isdiscretecache(integrator.cache) &&
        error("Derivatives are not defined for this stepper.")
    return if isfsal(integrator.alg) &&
            !has_stiff_interpolation(integrator.alg)
        # Special stiff interpolations do not store the
        # right value in fsallast
        integrator.fsallast
    else
        integrator(integrator.t, Val{1})
    end
end

@inline function SciMLBase.get_du!(out, integrator::ODEIntegratorType)
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

function u_modified!(integrator::ODEIntegratorType, bool::Bool)
    return integrator.u_modified = bool
end

function get_proposed_dt(integrator::ODEIntegratorType)
    return ifelse(integrator.opts.adaptive, integrator.dtpropose, integrator.dtcache)
end
function set_proposed_dt!(integrator::ODEIntegratorType, dt::Number)
    (integrator.dtpropose = dt; integrator.dtcache = dt)
end

function set_proposed_dt!(integrator::ODEIntegratorType, integrator2::ODEIntegratorType)
    integrator.dtpropose = integrator2.dtpropose
    integrator.dtcache = integrator2.dtcache
    integrator.qold = integrator2.qold
    integrator.erracc = integrator2.erracc
    return integrator.dtacc = integrator2.dtacc
end

#TODO: Bigger caches for most algorithms
@inline function SciMLBase.get_tmp_cache(integrator::ODEIntegratorType)
    return get_tmp_cache(integrator, integrator.alg, integrator.cache)
end

# the ordering of the cache arrays is important!!!
@inline function SciMLBase.get_tmp_cache(
        integrator, alg::OrdinaryDiffEqAlgorithm,
        cache::OrdinaryDiffEqConstantCache
    )
    return nothing
end
@inline function SciMLBase.get_tmp_cache(
        integrator, alg::OrdinaryDiffEqAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp,)
end
@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.nlsolver.tmp, cache.atmp)
end
@inline function SciMLBase.get_tmp_cache(
        integrator, alg::OrdinaryDiffEqNewtonAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.nlsolver.tmp, cache.nlsolver.z)
end
@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.linsolve_tmp)
end

@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.utilde)
end
@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.dz)
end
@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqLinearExponentialAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp,)
end
@inline function SciMLBase.get_tmp_cache(
        integrator, alg::CompositeAlgorithm,
        cache::CompositeCache
    )
    return get_tmp_cache(integrator, alg.algs[1], cache.caches[1])
end
@inline function SciMLBase.get_tmp_cache(
        integrator, alg::CompositeAlgorithm,
        cache::DefaultCacheType
    )
    init_ith_default_cache(cache, alg.algs, cache.current)
    return if cache.current == 1
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

@inline function SciMLBase.get_tmp_cache(
        integrator, alg::DAEAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.nlsolver.cache.dz, cache.atmp)
end

function full_cache(integrator::ODEIntegratorType)
    # for DefaultCache, we need to make sure to initialize all the caches in case they get switched to later
    if integrator.cache isa DefaultCacheType
        (; alg, cache) = integrator
        algs = alg.algs
        init_ith_default_cache(cache, algs, 1)
        init_ith_default_cache(cache, algs, 2)
        init_ith_default_cache(cache, algs, 3)
        init_ith_default_cache(cache, algs, 4)
        init_ith_default_cache(cache, algs, 5)
        init_ith_default_cache(cache, algs, 6)
    end
    return full_cache(integrator.cache)
end
function full_cache(cache::CompositeCache)
    return Iterators.flatten(full_cache(c) for c in cache.caches)
end
function full_cache(cache::DefaultCacheType)
    caches = (
        cache.cache1, cache.cache2, cache.cache3, cache.cache4, cache.cache5, cache.cache6,
    )
    return Iterators.flatten(full_cache(c) for c in caches)
end

function SciMLBase.add_tstop!(integrator::ODEIntegratorType, t)
    integrator.tdir * (t - integrator.t) < zero(integrator.t) &&
        error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
    return push!(integrator.opts.tstops, integrator.tdir * t)
end

SciMLBase.has_tstop(integrator::ODEIntegratorType) = !isempty(integrator.opts.tstops)
SciMLBase.first_tstop(integrator::ODEIntegratorType) = first(integrator.opts.tstops)
SciMLBase.pop_tstop!(integrator::ODEIntegratorType) = pop!(integrator.opts.tstops)

function SciMLBase.add_saveat!(integrator::ODEIntegratorType, t)
    integrator.tdir * (t - integrator.t) < zero(integrator.t) &&
        error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
    return push!(integrator.opts.saveat, integrator.tdir * t)
end

function resize!(integrator::ODEIntegratorType, i::Int)
    (; cache) = integrator

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
    return resize_non_user_cache!(integrator, cache, i)
end
# we can't use resize!(..., i::Union{Int, NTuple{N,Int}}) where {N} because of method ambiguities with DiffEqBase
function resize!(integrator::ODEIntegratorType, i::NTuple{N, Int}) where {N}
    (; cache) = integrator

    for c in full_cache(cache)
        resize!(c, i)
    end
    !isnothing(integrator.fsalfirst) && resize!(integrator.fsalfirst, i)
    !isnothing(integrator.fsallast) && resize!(integrator.fsallast, i)
    resize_f!(integrator.f, i)
    # TODO the parts below need to be adapted for implicit methods
    isdefined(integrator.cache, :nlsolver) && resize_nlsolver!(integrator, i)
    resize_J_W!(cache, integrator, i)
    return resize_non_user_cache!(integrator, cache, i)
end

# default fallback
resize_f!(f, i) = nothing

function resize_f!(f::SplitFunction, i)
    resize!(f._func_cache, i)
    return nothing
end

function resize_J_W! end

function resize_non_user_cache!(integrator::ODEIntegratorType, i::Int)
    return resize_non_user_cache!(integrator, integrator.cache, i)
end
function deleteat_non_user_cache!(integrator::ODEIntegratorType, i)
    return deleteat_non_user_cache!(integrator, integrator.cache, i)
end
function addat_non_user_cache!(integrator::ODEIntegratorType, i)
    return addat_non_user_cache!(integrator, integrator.cache, i)
end

resize_non_user_cache!(integrator::ODEIntegratorType, cache, i) = nothing

function resize_non_user_cache!(integrator::ODEIntegratorType, cache::CompositeCache, i)
    for _cache in cache.caches
        resize_non_user_cache!(integrator, _cache, i)
    end
    return
end

function deleteat_non_user_cache!(integrator::ODEIntegratorType, cache::CompositeCache, i)
    for _cache in cache.caches
        deleteat_non_user_cache!(integrator, _cache, i)
    end
    return
end

function addat_non_user_cache!(integrator::ODEIntegratorType, cache::CompositeCache, i)
    for _cache in cache.caches
        addat_non_user_cache!(integrator, _cache, i)
    end
    return
end

function deleteat_non_user_cache!(integrator::ODEIntegratorType, cache, idxs)
    # ordering doesn't matter in deterministic cache, so just resize
    # to match the size of u
    i = length(integrator.u)
    return resize_non_user_cache!(integrator, cache, i)
end

function addat_non_user_cache!(integrator::ODEIntegratorType, cache, idxs)
    # ordering doesn't matter in deterministic cache, so just resize
    # to match the size of u
    i = length(integrator.u)
    return resize_non_user_cache!(integrator, cache, i)
end

function deleteat!(integrator::ODEIntegratorType, idxs)
    for c in full_cache(integrator)
        deleteat!(c, idxs)
    end
    return deleteat_non_user_cache!(integrator, integrator.cache, idxs)
end

function addat!(integrator::ODEIntegratorType, idxs)
    for c in full_cache(integrator)
        addat!(c, idxs)
    end
    return addat_non_user_cache!(integrator, integrator.cache, idxs)
end

function terminate!(integrator::ODEIntegratorType, retcode = ReturnCode.Terminated)
    integrator.sol = SciMLBase.solution_new_retcode(integrator.sol, retcode)
    return integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end

const EMPTY_ARRAY_OF_PAIRS = Pair[]

SciMLBase.has_reinit(integrator::ODEIntegratorType) = true
function SciMLBase.reinit!(
        integrator::ODEIntegratorType, u0 = integrator.sol.prob.u0;
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
        reinit_retcode = true
    )
    if reinit_dae && SciMLBase.has_initializeprob(integrator.sol.prob.f)
        # This is `remake` infrastructure. `reinit!` is somewhat like `remake` for
        # integrators, so we reuse some of the same pieces. If we pass `integrator.p`
        # for `p`, it means we don't want to change it. If we pass `missing`, this
        # function may (correctly) assume `newp` aliases `prob.p` and copy it, which we
        # want to avoid. So we pass an empty array of pairs to make it think this is
        # a symbolic `remake` and it can modify `newp` inplace. The array of pairs is a
        # const global to avoid allocating every time this function is called.
        u0,
            newp = SciMLBase.late_binding_update_u0_p(
            integrator.sol.prob, u0,
            EMPTY_ARRAY_OF_PAIRS, t0, u0, integrator.p
        )
        if newp !== integrator.p
            integrator.p = newp
            sol = integrator.sol
            @reset sol.prob.p = newp
            integrator.sol = sol
        end
    end
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
    reinit_tstops!(tType, integrator.opts.tstops, tstops, d_discontinuities, tspan)
    reinit_saveat!(tType, integrator.opts.saveat, saveat, tspan)
    reinit_d_discontinuities!(tType, integrator.opts.d_discontinuities, d_discontinuities, tspan)
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

function SciMLBase.auto_dt_reset!(integrator::ODEIntegratorType)
    integrator.dt = ode_determine_initdt(
        integrator.u, integrator.t,
        integrator.tdir, integrator.opts.dtmax,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, integrator.sol.prob,
        integrator
    )
    integrator.dtpropose = integrator.dt
    return increment_nf!(integrator.stats, 2)
end

function increment_nf!(stats, amt = 1)
    return stats.nf += amt
end

function SciMLBase.set_t!(integrator::ODEIntegratorType, t::Real)
    if integrator.opts.save_everystep
        error(
            "Integrator time cannot be reset unless it is initialized",
            " with save_everystep=false"
        )
    end
    return if alg_extrapolates(integrator.alg) || !isdtchangeable(integrator.alg)
        reinit!(
            integrator, integrator.u;
            t0 = t,
            reset_dt = false,
            reinit_callbacks = false,
            reinit_cache = false
        )
    else
        integrator.t = t
    end
end

function SciMLBase.set_u!(integrator::ODEIntegratorType, u)
    if integrator.opts.save_everystep
        error(
            "Integrator state cannot be reset unless it is initialized",
            " with save_everystep=false"
        )
    end
    integrator.u = u
    return u_modified!(integrator, true)
end

SciMLBase.has_stats(i::ODEIntegratorType) = true

DiffEqBase.get_tstops(integ::ODEIntegratorType) = integ.opts.tstops
DiffEqBase.get_tstops_array(integ::ODEIntegratorType) = get_tstops(integ).valtree
DiffEqBase.get_tstops_max(integ::ODEIntegratorType) = maximum(get_tstops_array(integ))

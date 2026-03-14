@inline function DiffEqBase.change_t_via_interpolation!(
        integrator::Union{SDEIntegrator, AbstractSDDEIntegrator}, t,
        modify_save_endpoint::Type{Val{T}} = Val{false},
        reinitialize_alg = nothing
    ) where {T}
    # Can get rid of an allocation here with a function
    # get_tmp_arr(integrator.cache) which gives a pointer to some
    # cache array which can be modified.
    return if integrator.tdir * t < integrator.tdir * integrator.tprev
        error("Current interpolant only works between tprev and t")
    elseif t != integrator.t
        if integrator.u isa AbstractArray
            integrator(integrator.u, t)
        else
            integrator.u = integrator(t)
        end
        reject_step!(integrator, t - integrator.tprev) #this only changes dt and noise, so no interpolation problems
        integrator.t = t
        # dt should be the step from tprev to t (for correct interpolation in savevalues!)
        # This matches the behavior in OrdinaryDiffEq
        integrator.dt = integrator.t - integrator.tprev
        integrator.sqdt = sqrt(abs(integrator.dt))

        # reeval_internals_due_to_modification!(integrator) # Not necessary for linear interp
        if T
            solution_endpoint_match_cur_integrator!(integrator)
        end
    end
end

function (integrator::SDEIntegrator)(t, deriv::Type = Val{0}; idxs = nothing)
    return current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::SDEIntegrator)(
        val::AbstractArray, t::Union{Number, AbstractArray},
        deriv::Type = Val{0}; idxs = nothing
    )
    return current_interpolant!(val, t, integrator, idxs, deriv)
end

# set_proposed_dt!(::SDEIntegrator, ::SDEIntegrator) now provided by ODE's
# set_proposed_dt!(::ODEIntegrator, ::ODEIntegrator)

#TODO: Bigger caches for most algorithms
# avoid method ambiguity
for typ in (StochasticDiffEqAlgorithm, StochasticDiffEqNewtonAdaptiveAlgorithm)
    @eval @inline DiffEqBase.get_tmp_cache(
        integrator::SDEIntegrator, alg::$typ,
        cache::StochasticDiffEqConstantCache
    ) = nothing
end
@inline DiffEqBase.get_tmp_cache(integrator::SDEIntegrator, alg, cache) = (cache.tmp,)
@inline DiffEqBase.get_tmp_cache(
    integrator::SDEIntegrator, alg::StochasticDiffEqNewtonAdaptiveAlgorithm,
    cache
) = (cache.nlsolver.tmp, cache.nlsolver.ztmp)
@inline DiffEqBase.get_tmp_cache(
    integrator::SDEIntegrator, alg::StochasticCompositeAlgorithm,
    cache
) = get_tmp_cache(integrator, alg.algs[1], cache.caches[1])

function full_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(full_cache(c) for c in integrator.caches)
end

ratenoise_cache(integrator::SDEIntegrator) = ratenoise_cache(integrator.cache)
function ratenoise_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(ratenoise_cache(c) for c in integrator.caches)
end

rand_cache(integrator::SDEIntegrator) = rand_cache(integrator.cache)
function rand_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(rand_cache(c) for c in integrator.caches)
end

jac_iter(integrator::SDEIntegrator) = jac_iter(integrator.cache)
function jac_iter(integrator::StochasticCompositeCache)
    return Iterators.flatten(jac_iter(c) for c in integrator.caches)
end

resize!(integrator::SDEIntegrator, i::Int) = resize!(integrator, integrator.cache, i)

function resize!(integrator::SDEIntegrator, cache, i)
    # This has to go first!
    resize_non_user_cache!(integrator, cache, i)
    for c in full_cache(integrator)
        resize!(c, i)
    end
    for c in ratenoise_cache(integrator)
        resize!(c, i)
    end
    return
end

function resize_noise!(integrator, cache, bot_idx, i)
    for c in integrator.W.S₁
        resize!(c[2], i)
        if alg_needs_extra_process(integrator.alg)
            resize!(c[3], i)
        end
        if i >= bot_idx # fill in rands
            fill_new_noise_caches!(integrator, c, c[1], bot_idx:i)
        end
    end
    for c in integrator.W.S₂
        resize!(c[2], i)
        if alg_needs_extra_process(integrator.alg)
            resize!(c[3], i)
        end
        if i >= bot_idx # fill in rands
            fill_new_noise_caches!(integrator, c, c[1], bot_idx:i)
        end
    end
    resize!(integrator.W.dW, i)
    integrator.W.dW[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dWtilde, i)
    integrator.W.dWtilde[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dWtmp, i)
    integrator.W.dWtmp[end] = zero(eltype(integrator.u))
    resize!(integrator.W.curW, i)
    integrator.W.curW[end] = zero(eltype(integrator.u))
    DiffEqNoiseProcess.resize_stack!(integrator.W, i)

    if alg_needs_extra_process(integrator.alg)
        resize!(integrator.W.dZ, i)
        integrator.W.dZ[end] = zero(eltype(integrator.u))
        resize!(integrator.W.dZtilde, i)
        integrator.W.dZtilde[end] = zero(eltype(integrator.u))
        resize!(integrator.W.dZtmp, i)
        integrator.W.dZtmp[end] = zero(eltype(integrator.u))
        resize!(integrator.W.curZ, i)
        integrator.W.curZ[end] = zero(eltype(integrator.u))
    end
    return if i >= bot_idx # fill in rands
        fill!(@view(integrator.W.curW[bot_idx:i]), zero(eltype(integrator.u)))
        if alg_needs_extra_process(integrator.alg)
            fill!(@view(integrator.W.curZ[bot_idx:i]), zero(eltype(integrator.u)))
        end
    end
end

@inline function fill_new_noise_caches!(integrator, c, scaling_factor, idxs)
    return if isinplace(integrator.W)
        integrator.W.dist(
            @view(c[2][idxs]), integrator.W, scaling_factor,
            integrator.u, integrator.p, integrator.t, integrator.W.rng
        )
        if alg_needs_extra_process(integrator.alg)
            integrator.W.dist(
                @view(c[3][idxs]), integrator.W, scaling_factor,
                integrator.u, integrator.p, integrator.t, integrator.W.rng
            )
        end
    else
        c[2][idxs] .= integrator.noise(length(idxs), integrator, scaling_factor)
        if alg_needs_extra_process(integrator.alg)
            c[3][idxs] .= integrator.noise(length(idxs), integrator, scaling_factor)
        end
    end
end

function resize_non_user_cache!(integrator::SDEIntegrator, cache, i)
    bot_idx = length(integrator.u) + 1
    return if is_diagonal_noise(integrator.sol.prob)
        resize_noise!(integrator, cache, bot_idx, i)
        for c in rand_cache(integrator)
            resize!(c, i)
        end
    end
end

function deleteat!(integrator::SDEIntegrator, idxs)
    deleteat_non_user_cache!(integrator, integrator.cache, idxs)
    for c in full_cache(integrator)
        deleteat!(c, idxs)
    end
    for c in ratenoise_cache(integrator)
        deleteat!(c, idxs)
    end
    return
end

function addat!(integrator::SDEIntegrator, idxs)
    addat_non_user_cache!(integrator, integrator.cache, idxs)
    for c in full_cache(integrator)
        addat!(c, idxs)
    end
    for c in ratenoise_cache(integrator)
        addat!(c, idxs)
    end
    return
end

function deleteat_non_user_cache!(integrator::SDEIntegrator, cache, idxs)
    return if is_diagonal_noise(integrator.sol.prob)
        deleteat_noise!(integrator, cache, idxs)
        for c in rand_cache(integrator)
            deleteat!(c, idxs)
        end
    end
end

function addat_non_user_cache!(integrator::SDEIntegrator, cache, idxs)
    return if is_diagonal_noise(integrator.sol.prob)
        addat_noise!(integrator, cache, idxs)
        for c in rand_cache(integrator)
            addat!(c, idxs)
        end
    end
end

function deleteat_noise!(integrator, cache, idxs)
    for c in integrator.W.S₁
        deleteat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            deleteat!(c[3], idxs)
        end
    end
    for c in integrator.W.S₂
        deleteat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            deleteat!(c[3], idxs)
        end
    end
    deleteat!(integrator.W.dW, idxs)
    deleteat!(integrator.W.dWtilde, idxs)
    deleteat!(integrator.W.dWtmp, idxs)
    deleteat!(integrator.W.curW, idxs)
    DiffEqNoiseProcess.resize_stack!(integrator.W, length(integrator.u))

    return if alg_needs_extra_process(integrator.alg)
        deleteat!(integrator.W.curZ, idxs)
        deleteat!(integrator.W.dZtmp, idxs)
        deleteat!(integrator.W.dZtilde, idxs)
        deleteat!(integrator.W.dZ, idxs)
    end
end

function addat_noise!(integrator, cache, idxs)
    for c in integrator.W.S₁
        addat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            addat!(c[3], idxs)
        end
        fill_new_noise_caches!(integrator, c, c[1], idxs)
    end
    for c in integrator.W.S₂
        addat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            addat!(c[3], idxs)
        end
        fill_new_noise_caches!(integrator, c, c[1], idxs)
    end

    addat!(integrator.W.dW, idxs)
    integrator.W.dW[idxs] .= zero(eltype(integrator.u))
    addat!(integrator.W.curW, idxs)
    integrator.W.curW[idxs] .= zero(eltype(integrator.u))
    if alg_needs_extra_process(integrator.alg)
        addat!(integrator.W.dZ, idxs)
        integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
        addat!(integrator.W.curZ, idxs)
        integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
    end

    i = length(integrator.u)
    resize!(integrator.W.dWtilde, i)
    resize!(integrator.W.dWtmp, i)
    DiffEqNoiseProcess.resize_stack!(integrator.W, i)
    if alg_needs_extra_process(integrator.alg)
        resize!(integrator.W.dZtmp, i)
        resize!(integrator.W.dZtilde, i)
    end

    # fill in rands
    fill!(@view(integrator.W.curW[idxs]), zero(eltype(integrator.u)))
    return if alg_needs_extra_process(integrator.alg)
        fill!(@view(integrator.W.curZ[idxs]), zero(eltype(integrator.u)))
    end
end

function DiffEqBase.reinit!(
        integrator::SDEIntegrator, u0 = integrator.sol.prob.u0;
        t0 = integrator.sol.prob.tspan[1], tf = integrator.sol.prob.tspan[2],
        erase_sol = true,
        tstops = integrator.opts.tstops_cache,
        saveat = integrator.opts.saveat_cache,
        d_discontinuities = integrator.opts.d_discontinuities_cache,
        reinit_cache = true, reinit_callbacks = true,
        initialize_save = true,
        reset_dt = (integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive,
        rng = nothing
    )
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.u, u0)
        recursivecopy!(integrator.uprev, integrator.u)
    else
        integrator.u = u0
        integrator.uprev = integrator.u
    end

    integrator.t = t0
    integrator.tprev = t0

    tType = typeof(integrator.t)
    tspan = (tType(t0), tType(tf))

    # Stash callable tstops (e.g. SymbolicTstops) and use empty tuple for heap init.
    if tstops isa AbstractArray || tstops isa Tuple || tstops isa Number
        _tstops_callable = nothing
    else
        _tstops_callable = tstops
        tstops = ()
    end

    integrator.opts.tstops = OrdinaryDiffEqCore.initialize_tstops(tType, tstops, d_discontinuities, tspan)
    integrator.opts.saveat = OrdinaryDiffEqCore.initialize_saveat(tType, saveat, tspan)
    integrator.opts.d_discontinuities = OrdinaryDiffEqCore.initialize_d_discontinuities(
        tType, d_discontinuities, tspan
    )

    if erase_sol
        if integrator.opts.save_start
            resize_start = 1
        else
            resize_start = 0
        end
        resize!(integrator.sol.u, resize_start)
        resize!(integrator.sol.t, resize_start)
        if integrator.sol.u_analytic !== nothing
            resize!(integrator.sol.u_analytic, 0)
        end
        if integrator.alg isa StochasticDiffEqCompositeAlgorithm
            resize!(integrator.sol.alg_choice, resize_start)
        end
        integrator.saveiter = resize_start
    end
    integrator.iter = 0
    integrator.success_iter = 0

    # full re-initialize the PI in timestepping
    integrator.qold = integrator.opts.qoldinit
    integrator.q11 = typeof(integrator.t)(1)

    if rng !== nothing
        SciMLBase.set_rng!(integrator, rng)
    end

    if reset_dt
        auto_dt_reset!(integrator)
    end

    if reinit_callbacks
        initialize_callbacks!(integrator, initialize_save)
    end

    if reinit_cache
        initialize!(integrator, integrator.cache)
    end

    # Evaluate callable tstops now that callbacks are re-initialized.
    if _tstops_callable !== nothing
        for ts in _tstops_callable(integrator.p, tspan)
            add_tstop!(integrator, ts)
        end
    end

    return reinit!(integrator.W, integrator.dt)
end

function DiffEqBase.auto_dt_reset!(integrator::SDEIntegrator)
    return integrator.dt = sde_determine_initdt(
        integrator.u, integrator.t,
        integrator.tdir, integrator.opts.dtmax, integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, integrator.sol.prob, get_current_alg_order(
            integrator.alg, integrator.cache
        ),
        integrator
    )
end

# get_du and get_du! now provided by ODE's SciMLBase.get_du(::ODEIntegrator).
# For SDE (isfsal=false), ODE's version falls through to integrator(t, Val{1})
# which calls SDE's linear interpolation, returning (u - uprev) / dt.

function DiffEqBase.set_t!(integrator::SDEIntegrator, t::Real)
    if integrator.opts.save_everystep
        error(
            "Integrator time cannot be reset unless it is initialized",
            " with save_everystep=false"
        )
    end
    return if !isdtchangeable(integrator.alg)
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

"""
    SciMLBase.set_rng!(integrator::SDEIntegrator, rng) -> nothing

Replace the integrator's random number generator. The new RNG must be the same
concrete type as the current one (the type is baked into the integrator's type
parameters).

## Noise process behavior

- **Framework-constructed noise** (no `noise` on the problem): `set_rng!` also
  updates `W.rng` and `P.rng` so that all framework-managed randomness uses the
  new RNG. After the call, `integrator.rng === integrator.W.rng`.
- **User-provided noise** (`SDEProblem(...; noise = my_W)`): `set_rng!` only
  updates `integrator.rng`. The noise process keeps its own RNG — the user is
  responsible for managing it (e.g., via `Random.seed!(integrator.W.rng, seed)`
  or replacing it directly with `integrator.W.rng = new_rng`).

## See also

[`reinit!`](@ref) accepts an `rng` keyword that delegates to `set_rng!`.
"""
function SciMLBase.set_rng!(integrator::SDEIntegrator, rng)
    R = typeof(integrator.rng)
    if !isa(rng, R)
        throw(
            ArgumentError(
                "Cannot set RNG of type $(typeof(rng)) on an integrator " *
                    "whose RNG type parameter is $R. " *
                    "Construct a new integrator via `init(prob, alg; rng = your_rng)` instead."
            )
        )
    end
    integrator.rng = rng
    # Sync framework-constructed noise processes only
    if integrator.noise === nothing && integrator.W !== nothing
        integrator.W.rng = rng
    end
    # P (CompoundPoissonProcess) is always framework-constructed when present
    if integrator.P !== nothing
        integrator.P.rng = rng
    end
    return nothing
end

function _get_alias_noise_from_kwargs(; alias_noise = nothing, alias = nothing, kwargs...)
    if alias_noise !== nothing
        return alias_noise
    elseif alias !== nothing && hasproperty(alias, :alias_noise) && alias.alias_noise !== nothing
        return alias.alias_noise
    else
        return true
    end
end

function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractRODEProblem,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        kwargs...
    )
    integrator = DiffEqBase.__init(prob, alg; kwargs...)
    solve!(integrator)
    if prob isa DiffEqBase.AbstractRODEProblem &&
            typeof(prob.noise) == typeof(integrator.sol.W) &&
            _get_alias_noise_from_kwargs(; kwargs...)
        copy!(prob.noise, integrator.sol.W)
    end
    return integrator.sol
end

# More specific method for JumpProblem to win over JumpProcesses.jl's ambiguity fix dispatch
function DiffEqBase.__solve(
        prob::JumpProblem,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        merge_callbacks = true, kwargs...
    )
    kwargs = DiffEqBase.merge_problem_kwargs(prob; merge_callbacks, kwargs...)
    integrator = DiffEqBase.__init(prob, alg; kwargs...)
    solve!(integrator)
    if concrete_prob(prob) isa DiffEqBase.AbstractRODEProblem &&
            typeof(concrete_prob(prob).noise) == typeof(integrator.sol.W) &&
            _get_alias_noise_from_kwargs(; kwargs...)
        copy!(concrete_prob(prob).noise, integrator.sol.W)
    end
    return integrator.sol
end

# Make it easy to grab the RODEProblem/SDEProblem/DiscreteProblem from the keyword arguments
concrete_prob(prob) = prob
concrete_prob(prob::JumpProblem) = prob.prob

"""
    _resolve_rng(rng, seed, prob) -> (rng, seed, rng_provided)

Resolve the RNG and seed for an SDE/RODE integration from the user-provided
`rng` and `seed` kwargs plus the problem's stored seed.

## Priority

1. **`rng` provided**: Use it directly. `rng_provided = true`. If the RNG is a
   `TaskLocalRNG`, it is converted to a concrete `Xoshiro` seeded from a draw
   from the task-local stream, and the derived seed is returned as `seed`. This
   avoids type mismatches (DiffEqNoiseProcess copies `TaskLocalRNG` into `Xoshiro`
   internally) and prevents duplicate random streams. For other RNG types,
   `seed = UInt64(0)` (sentinel for "no seed").
2. **`seed` provided (nonzero)**: Construct `Xoshiro(seed)`. `rng_provided = false`.
3. **Problem has a seed** (`prob.seed != 0`): Use it. `rng_provided = false`.
4. **Neither**: Generate a random seed and construct `Xoshiro` from it.
   `rng_provided = false`.

## Return values

- `rng::AbstractRNG` — the RNG to use for noise processes and `integrator.rng`.
- `seed::UInt64` — the seed for solution metadata (`RODESolution.seed::UInt64`).
  For `TaskLocalRNG` conversion, this is the derived seed. For other user-provided
  RNGs, `UInt64(0)` (no seed was used).
- `rng_provided::Bool` — `true` when the user explicitly passed an `rng` kwarg.
  Controls whether downstream consumers (JumpProcesses reseeding, noise process
  reseeding) should skip their own seed-based initialization.
"""
function _resolve_rng(rng, seed, prob)
    if rng !== nothing
        if !(rng isa Random.AbstractRNG)
            throw(
                ArgumentError(
                    "`rng` must be an `AbstractRNG`, got $(typeof(rng))."
                )
            )
        end
        # TaskLocalRNG is a zero-field singleton pointing to shared task-local
        # state. Storing it directly would cause type mismatches with noise
        # processes (DiffEqNoiseProcess converts it to Xoshiro internally) and
        # means set_rng! cannot sync integrator.rng with W.rng. Convert to a
        # concrete Xoshiro seeded from one draw of the task-local stream.
        if rng isa Random.TaskLocalRNG
            _seed = rand(rng, UInt64)
            return Random.Xoshiro(_seed), _seed, true
        end
        return rng, UInt64(0), true
    end
    _seed = if iszero(seed)
        if (!(prob isa DiffEqBase.AbstractRODEProblem) || iszero(prob.seed))
            seed_multiplier() * rand(UInt64)
        else
            prob.seed
        end
    else
        seed
    end
    return Random.Xoshiro(_seed), _seed, false
end

# Quick fix: replace ::Union{DiffEqBase.AbstractRODEProblem, JumpProblem} with a
# concrete arguments using an eval loop to avoid ambiguities.
for Prob ∈ (DiffEqBase.AbstractRODEProblem, JumpProblem)
    @eval begin
# rng kwarg: Pre-constructed AbstractRNG for all framework-managed randomness
# (noise processes, integrator.rng). Takes priority over `seed` when provided.
# When omitted, an Xoshiro is constructed from `seed` (or a random seed if
# `seed` is also omitted). Only controls framework-constructed randomness;
# user-provided noise processes (`prob.noise`) keep their own RNG.
function DiffEqBase.__init(
        _prob::$Prob,
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        saveat = (),
        tstops = (),
        d_discontinuities = (),
        save_idxs = nothing,
        save_everystep = isempty(saveat),
        save_noise = save_everystep && (
            typeof(concrete_prob(_prob).f) <: Tuple ?
                DiffEqBase.has_analytic(concrete_prob(_prob).f[1]) :
                DiffEqBase.has_analytic(concrete_prob(_prob).f)
        ),
        save_on = true,
        save_start = save_everystep || isempty(saveat) || saveat isa Number ? true :
            concrete_prob(_prob).tspan[1] in saveat,
        save_end = nothing,
        callback = nothing,
        dense = false, # save_everystep && isempty(saveat),
        calck = false, #(!isempty(setdiff(saveat,tstops)) || dense),
        dt = eltype(concrete_prob(_prob).tspan)(0),
        adaptive = isadaptive(_prob, alg),
        gamma = isadaptive(alg) ? 9 // 10 : 0,
        abstol = nothing,
        reltol = nothing,
        qmin = qmin_default(alg),
        qmax = qmax_default(alg),
        qsteady_min = qsteady_min_default(alg),
        qsteady_max = qsteady_max_default(alg),
        beta2 = nothing,
        beta1 = nothing,
        qoldinit = isadaptive(alg) ? 1 // 10^4 : 0,
        controller = nothing,
        fullnormalize = true,
        failfactor = 2,
        delta = delta_default(alg),
        maxiters = adaptive ? 1000000 : typemax(Int),
        dtmax = eltype(concrete_prob(_prob).tspan)((concrete_prob(_prob).tspan[end] - concrete_prob(_prob).tspan[1])),
        dtmin = DiffEqBase.prob2dtmin(concrete_prob(_prob)),
        internalnorm = ODE_DEFAULT_NORM,
        isoutofdomain = ODE_DEFAULT_ISOUTOFDOMAIN,
        unstable_check = ODE_DEFAULT_UNSTABLE_CHECK,
        verbose = true, force_dtmin = false,
        timeseries_errors = true, dense_errors = false,
        advance_to_tstop = false, stop_at_next_tstop = false,
        initialize_save = true,
        progress = false, progress_steps = 1000, progress_name = "SDE",
        progress_message = ODE_DEFAULT_PROG_MESSAGE,
        progress_id = :StochasticDiffEq,
        userdata = nothing,
        save_discretes = true,
        initialize_integrator = true,
        seed = UInt64(0),
        rng = nothing,
        alias = nothing,
        initializealg = OrdinaryDiffEqCore.DefaultInit(),
        kwargs...
    )
    # Merge JumpProblem kwargs (tstops, callbacks, etc.) for direct init() callers.
    if _prob isa JumpProblem
        kwargs = DiffEqBase.merge_problem_kwargs(_prob; kwargs...)
    end

    is_sde = _prob isa SDEProblem

    use_old_kwargs = haskey(kwargs, :alias_u0) || haskey(kwargs, :alias_jumps) ||
        haskey(kwargs, :alias_noise)

    if use_old_kwargs
        aliases = ODEAliasSpecifier()
        if haskey(kwargs, :alias_u0)
            message = "`alias_u0` keyword argument is deprecated, to set `alias_u0`,
            please use an SDEAliasSpecifier or RODEAliasSpecifier, e.g. `solve(prob, alias = SDEAliasSpecifier(alias_u0 = true))`"
            Base.depwarn(message, :init)
            Base.depwarn(message, :solve)
            alias_u0 = values(kwargs).alias_u0
        else
            alias_u0 = nothing
        end

        if haskey(kwargs, :alias_jumps)
            message = "`alias_jumps` keyword argument is deprecated, to set `alias_jumps`,
            please use an SDEAliasSpecifier or RODEAliasSpecifier, e.g. `solve(prob, alias = SDEAliasSpecifier(alias_jumps = true))`"
            Base.depwarn(message, :init)
            Base.depwarn(message, :solve)
            alias_jumps = values(kwargs).alias_jumps
        else
            alias_jumps = nothing
        end

        if haskey(kwargs, :alias_noise)
            message = "`alias_noise` keyword argument is deprecated, to set `alias_noise`,
            please use a RODEAliasSpecifier, e.g. `solve(prob, alias = RODEAliasSpecifier(alias_noise = true))`"
            Base.depwarn(message, :init)
            Base.depwarn(message, :solve)
            alias_noise = values(kwargs).alias_noise
        else
            alias_noise = nothing
        end

        aliases = if alias_noise !== nothing
            SciMLBase.RODEAliasSpecifier(; alias_u0, alias_jumps, alias_noise)
        elseif is_sde
            SciMLBase.SDEAliasSpecifier(; alias_u0, alias_jumps)
        else
            SciMLBase.RODEAliasSpecifier(; alias_u0, alias_jumps)
        end

    else
        # If alias isa Bool, all fields of SDEAliasSpecifier set to alias
        if alias isa Bool
            aliases = is_sde ? SciMLBase.SDEAliasSpecifier(; alias) :
                SciMLBase.RODEAliasSpecifier(; alias)
        elseif alias isa SciMLBase.SDEAliasSpecifier
            aliases = alias
        elseif alias isa SciMLBase.RODEAliasSpecifier
            aliases = alias
        else
            # Default case (including isnothing(alias))
            aliases = is_sde ? SciMLBase.SDEAliasSpecifier() :
                SciMLBase.RODEAliasSpecifier()
        end
    end

    prob = concrete_prob(_prob)

    _rng, _seed, _rng_provided = _resolve_rng(rng, seed, prob)

    if _prob isa JumpProblem
        alias_jumps = isnothing(aliases.alias_jumps) ? Threads.threadid() == 1 :
            aliases.alias_jumps
        _jump_seed = _rng_provided ? nothing : _seed
        if !alias_jumps
            _prob = JumpProcesses.resetted_jump_problem(_prob, _jump_seed)
        else
            JumpProcesses.reset_jump_problem!(_prob, _jump_seed)
        end
    end

    # Grab the deepcopied version for caching reasons.
    prob = concrete_prob(_prob)

    if typeof(prob.f) <: Tuple
        if any(mm != I for mm in prob.f.mass_matrix)
            error("This solver is not able to use mass matrices.")
        end
    elseif prob isa DiffEqBase.AbstractRODEProblem && prob.f.mass_matrix != I &&
            !alg_mass_matrix_compatible(alg)
        error("This solver is not able to use mass matrices.")
    end

    if !isempty(saveat) && dense
        @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
    end

    if prob isa DiffEqBase.AbstractRODEProblem && typeof(prob.noise) <: NoiseProcess &&
            prob.noise.bridge === nothing && adaptive
        error("Bridge function must be given for adaptivity. Either declare this function in noise process or set adaptive=false")
    end

    if !alg_compatible(_prob, alg)
        error("The algorithm is not compatible with the chosen noise type. Please see the documentation on the solver methods")
    end

    if adaptive && !isadaptive(_prob, alg)
        error("The given solver is a Fixed timestep method and does not support adaptivity.")
    end

    progress && @logmsg(LogLevel(-1), progress_name, _id = progress_id, progress = 0)

    tType = eltype(prob.tspan)
    noise = prob isa DiffEqBase.AbstractRODEProblem ? prob.noise : nothing
    tspan = prob.tspan
    tdir = sign(tspan[end] - tspan[1])

    t = tspan[1]

    if !adaptive && iszero(dt) &&
            (tstops isa AbstractArray || tstops isa Tuple) && isempty(tstops)
        error("Fixed timestep methods require a choice of dt or choosing the tstops")
    end

    if alg isa StochasticCompositeAlgorithm && alg.choice_function isa AutoSwitch
        auto = alg.choice_function
        alg = StochasticCompositeAlgorithm(
            alg.algs,
            AutoSwitchCache(
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

    g = prob isa DiffEqBase.AbstractSDEProblem ? prob.g : nothing

    if prob.u0 isa Tuple
        u = ArrayPartition(prob.u0, Val{true})
    end

    uType = typeof(u)
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    ks = ()

    uEltypeNoUnits = recursive_unitless_eltype(u)
    tTypeNoUnits = typeof(one(tType))

    if abstol === nothing
        if uBottomEltypeNoUnits === uBottomEltype && !(uBottomEltype <: Integer)
            abstol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        elseif uBottomEltype <: Integer
            abstol_internal = real(oneunit(uBottomEltype) * 1 // 10^2)
        else
            abstol_internal = 0
        end
    else
        abstol_internal = real.(abstol)
    end

    if reltol === nothing
        if uBottomEltypeNoUnits === uBottomEltype && !(uBottomEltype <: Integer)
            reltol_internal = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        elseif uBottomEltype <: Integer
            reltol_internal = real(oneunit(uBottomEltype) * 1 // 10^2)
        else
            reltol_internal = 0
        end
    else
        reltol_internal = real.(reltol)
    end

    dtmax > zero(dtmax) && tdir < 0 && (dtmax *= tdir) # Allow positive dtmax, but auto-convert
    # dtmin is all abs => does not care about sign already.

    if isinplace(prob) && u isa AbstractArray && eltype(u) <: Number &&
            uBottomEltypeNoUnits == uBottomEltype # Could this be more efficient for other arrays?
        if !(u isa ArrayPartition)
            rate_prototype = recursivecopy(u)
        else
            rate_prototype = similar(
                u, typeof.(
                    oneunit.(recursive_bottom_eltype.(u.x)) ./
                        oneunit(tType)
                )...
            )
        end
    else
        if uBottomEltypeNoUnits == uBottomEltype
            rate_prototype = u
        else # has units!
            rate_prototype = u / oneunit(tType)
        end
    end
    rateType = typeof(rate_prototype) ## Can be different if united

    if prob.f isa DynamicalSDEFunction
        noise_rate_prototype = rate_prototype.x[1]
    elseif is_diagonal_noise(prob)
        noise_rate_prototype = rate_prototype
    elseif prob isa DiffEqBase.AbstractRODEProblem
        if prob isa DiffEqBase.AbstractSDEProblem
            noise_rate_prototype = copy(prob.noise_rate_prototype)
        else
            noise_rate_prototype = copy(prob.rand_prototype)
        end
    else
        noise_rate_prototype = nothing
    end

    # Stash callable tstops (e.g. SymbolicTstops) and use empty tuple for heap init.
    if tstops isa AbstractArray || tstops isa Tuple || tstops isa Number
        _tstops_callable = nothing
    else
        _tstops_callable = tstops
        tstops = ()
    end

    tstops_internal = OrdinaryDiffEqCore.initialize_tstops(tType, tstops, d_discontinuities, tspan)
    saveat_internal = OrdinaryDiffEqCore.initialize_saveat(tType, saveat, tspan)
    d_discontinuities_internal = OrdinaryDiffEqCore.initialize_d_discontinuities(
        tType, d_discontinuities, tspan
    )

    ### Algorithm-specific defaults ###
    save_idxs,
        saved_subsystem = SciMLBase.get_save_idxs_and_saved_subsystem(prob, save_idxs)
    # if save_idxs === nothing
    #   ksEltype = Vector{rateType}
    # else
    #   ks_prototype = rate_prototype[save_idxs]
    #   ksEltype = Vector{typeof(ks_prototype)}
    # end

    # Initialize timeseries and ts vectors
    u_initial = save_idxs === nothing ? u : u[save_idxs]
    if save_idxs === nothing
        timeseries = Vector{uType}()
    else
        timeseries = Vector{typeof(u_initial)}()
    end
    ts = Vector{tType}()

    alg_choice = alg isa StochasticDiffEqCompositeAlgorithm ? Int[] : ()

    if !adaptive && save_everystep && tspan[2] - tspan[1] != Inf
        iszero(dt) ? steps = length(tstops) :
            steps = ceil(Int, internalnorm((tspan[2] - tspan[1]) / dt, tspan[1]))
        sizehint!(timeseries, steps + 1)
        sizehint!(ts, steps + 1)
    elseif save_everystep
        sizehint!(timeseries, 50)
        sizehint!(ts, 50)
    elseif !isempty(saveat_internal)
        sizehint!(timeseries, length(saveat_internal) + 1)
        sizehint!(ts, length(saveat_internal) + 1)
    else
        sizehint!(timeseries, 2)
        sizehint!(ts, 2)
    end

    if save_start
        saveiter = 1 # Starts at 1 so first save is at 2
        copyat_or_push!(ts, 1, t)
        if save_idxs === nothing
            copyat_or_push!(timeseries, 1, u)
        else
            copyat_or_push!(timeseries, 1, u_initial, Val{false})
        end
        if alg isa StochasticDiffEqCompositeAlgorithm
            copyat_or_push!(alg_choice, 1, 1)
        end
    else
        saveiter = 0
    end

    QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits

    uprev = recursivecopy(u)

    if !(uType <: AbstractArray)
        rand_prototype = zero(u ./ u) # Strip units and type info
        randType = typeof(rand_prototype)
    else
        randElType = uBottomEltypeNoUnits # Strip units and type info
        if prob.f isa DynamicalSDEFunction
            rand_prototype = copy(noise_rate_prototype)
        elseif is_diagonal_noise(prob)
            if u isa SArray
                rand_prototype = zero(u) # TODO: Array{randElType} for units
            else
                rand_prototype = (u .- u) ./ sqrt(oneunit(t))
            end
        elseif prob isa DiffEqBase.AbstractSDEProblem
            if issparse(u)
                rand_prototype = adapt(
                    DiffEqBase.parameterless_type(u), zeros(randElType, size(noise_rate_prototype, 2))
                )
            else
                rand_prototype = false .* noise_rate_prototype[1, :]
            end
        elseif prob isa DiffEqBase.AbstractRODEProblem
            rand_prototype = copy(prob.rand_prototype)
        else
            rand_prototype = nothing
        end
        randType = typeof(rand_prototype) # Strip units and type info
    end

    if _prob isa JumpProblem
        callbacks_internal = CallbackSet(callback, _prob.jump_callback)
    else
        callbacks_internal = CallbackSet(callback)
    end

    max_len_cb = DiffEqBase.max_vector_callback_length(callbacks_internal)
    if max_len_cb isa DiffEqBase.VectorContinuousCallback
        callback_cache = DiffEqBase.CallbackCache(max_len_cb.len, uBottomEltype, uBottomEltype)
    else
        callback_cache = nothing
    end

    if prob isa DiffEqBase.AbstractRODEProblem && prob.noise === nothing
        rswm = isadaptive(alg) ? RSWM(adaptivealg = :RSwM3) : RSWM(adaptivealg = :RSwM1)
        if isinplace(prob)
            #if isadaptive(alg) || callback !== nothing
            if alg_needs_extra_process(alg)
                if alg === PL1WM()
                    m = length(rand_prototype)
                    rand_prototype2 = similar(rand_prototype, Int(m * (m - 1) / 2))
                    rand_prototype2 .= false
                    W = WienerProcess!(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                elseif alg isa RKMilGeneral
                    m = length(rand_prototype)
                    if alg.p != nothing
                        rand_prototype2 = similar(rand_prototype, Int(m + alg.p * m * 2))
                        rand_prototype2 .= false
                    else
                        rand_prototype2 = nothing
                    end
                    W = WienerProcess!(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                elseif alg isa W2Ito1
                    m = 2
                    rand_prototype2 = zeros(eltype(rand_prototype), Int(m))
                    W = WienerProcess!(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                else
                    W = WienerProcess!(
                        t, rand_prototype, rand_prototype,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                end
            else
                W = WienerProcess!(
                    t, rand_prototype,
                    save_everystep = save_noise,
                    rng = _rng
                )
            end
            #=
            else
              if alg_needs_extra_process(alg)
                W = SimpleWienerProcess!(t,rand_prototype,rand_prototype,
                                   save_everystep=save_noise,
                                   rng = _rng)
              else
                W = SimpleWienerProcess!(t,rand_prototype,
                                   save_everystep=save_noise,
                                   rng = _rng)
              end
            end
            =#
        else
            #if isadaptive(alg) || callback !== nothing
            if alg_needs_extra_process(alg)
                if alg === PL1WM()
                    m = length(rand_prototype)
                    if rand_prototype isa Number
                        rand_prototype2 = nothing
                    else
                        rand_prototype2 = similar(rand_prototype, Int(m * (m - 1) / 2))
                        rand_prototype2 .= false
                    end
                    W = WienerProcess(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                elseif alg isa RKMilGeneral
                    m = length(rand_prototype)
                    if rand_prototype isa Number || alg.p === nothing
                        rand_prototype2 = nothing
                    else
                        rand_prototype2 = similar(rand_prototype, Int(m + alg.p * m * 2))
                        rand_prototype2 .= false
                    end
                    W = WienerProcess(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                elseif alg isa W2Ito1
                    m = 2
                    rand_prototype2 = zeros(eltype(rand_prototype), Int(m))
                    W = WienerProcess(
                        t, rand_prototype, rand_prototype2,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                else
                    W = WienerProcess(
                        t, rand_prototype, rand_prototype,
                        save_everystep = save_noise,
                        rng = _rng
                    )
                end
            else
                W = WienerProcess(
                    t, rand_prototype,
                    save_everystep = save_noise,
                    rng = _rng
                )
            end
            #=
            else
              if alg_needs_extra_process(alg)
                W = SimpleWienerProcess(t,rand_prototype,rand_prototype,
                                   save_everystep=save_noise,
                                   rng = _rng)
              else
                W = SimpleWienerProcess(t,rand_prototype,
                                   save_everystep=save_noise,
                                   rng = _rng)
              end
            end
            =#
        end
    elseif prob isa DiffEqBase.AbstractRODEProblem
        _alias_noise = if hasproperty(aliases, :alias_noise) && aliases.alias_noise !== nothing
            aliases.alias_noise
        else
            true
        end
        W = _alias_noise ? copy(prob.noise) : prob.noise

        if alg_needs_extra_process(alg) && (!hasproperty(W, :dZ) || W.dZ === nothing)
            error("Higher order solver requires extra Brownian process Z. Thus `WienerProcess(t, W0)` is insufficient, you must use `WienerProcess(t, W0, Z0)` where `Z` is another Brownian process")
        end

        if W.reset
            # Reseed (skip when user explicitly provided an RNG — reseeding
            # would mutate the shared _rng object, affecting the integrator
            # and any other noise processes using it)
            if !_rng_provided && W isa Union{NoiseProcess, NoiseTransport} && W.reseed
                Random.seed!(W.rng, _seed)
            end
            if W.curt != t
                reinit!(W, t, t0 = t)
            end

        elseif W.curt != t
            error("Starting time in the noise process is not the starting time of the simulation. The noise process should be re-initialized for repeated use")
        end
    else # Only a jump problem
        @assert _prob isa JumpProblem
        W = nothing
    end

    if _prob isa JumpProblem && _prob.regular_jump !== nothing
        if !isnothing(_prob.regular_jump.mark_dist) == nothing # https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/250
            error("Mark distributions are currently not supported in SimpleTauLeaping")
        end

        jump_prototype = zeros(_prob.regular_jump.numjumps)
        c = _prob.regular_jump.c

        if isinplace(_prob.regular_jump)
            rate_constants = zeros(_prob.regular_jump.numjumps)
            _prob.regular_jump.rate(rate_constants, u ./ u, prob.p, tspan[1])
            P = CompoundPoissonProcess!(
                _prob.regular_jump.rate, t, jump_prototype,
                computerates = !alg_control_rate(alg) || !adaptive,
                save_everystep = save_noise,
                rng = _rng
            )
            alg_control_rate(alg) && adaptive &&
                P.cache.rate(P.cache.currate, u, p, tspan[1])
        else
            rate_constants = _prob.regular_jump.rate(u ./ u, prob.p, tspan[1])
            P = CompoundPoissonProcess(
                _prob.regular_jump.rate, t, jump_prototype,
                save_everystep = save_noise,
                computerates = !alg_control_rate(alg) || !adaptive,
                rng = _rng
            )
            alg_control_rate(alg) && adaptive &&
                (P.cache.currate = P.cache.rate(u, p, tspan[1]))
        end

    else
        jump_prototype = nothing
        c = nothing
        P = nothing
        rate_constants = nothing
    end

    dW, dZ = isnothing(W) ? (nothing, nothing) : (W.dW, W.dZ)

    # Convert verbose argument to DEVerbosity
    verbose_internal = if verbose isa Bool
        verbose ? DEVerbosity(Standard()) : DEVerbosity(None())
    elseif verbose isa AbstractVerbosityPreset
        DEVerbosity(verbose)
    elseif verbose isa DEVerbosity
        verbose
    else
        throw(ArgumentError("verbose must be a Bool, AbstractVerbosityPreset, or DEVerbosity"))
    end

    cache = alg_cache(
        alg, prob, u, dW, dZ, p, rate_prototype, noise_rate_prototype,
        jump_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, f, t, dt, Val{isinplace(_prob)}, verbose_internal
    )

    if _prob isa JumpProblem && prob isa DiscreteProblem && prob isa Integer
        id = DiffEqBase.ConstantInterpolation(ts, timeseries)
    else
        id = LinearInterpolationData(timeseries, ts)
    end

    save_end_user = save_end
    save_end = save_end === nothing ?
        save_everystep || isempty(saveat) || saveat isa Number ||
        prob.tspan[2] in saveat : save_end

    # Setting up the step size controller
    if (beta1 !== nothing || beta2 !== nothing) && controller !== nothing
        throw(
            ArgumentError(
                "Setting both the legacy PID parameters `beta1, beta2 = $((beta1, beta2))` and the `controller = $controller` is not allowed."
            )
        )
    end

    if (beta1 !== nothing || beta2 !== nothing)
        message = "Providing the legacy PID parameters `beta1, beta2` is deprecated. Use the keyword argument `controller` instead."
        Base.depwarn(message, :init)
        Base.depwarn(message, :solve)
    end

    if controller === nothing
        controller = default_controller(alg, cache, QT(qoldinit), beta1, beta2)
    end

    # Restore callable for tstops_cache so reinit! can re-evaluate it.
    if _tstops_callable !== nothing
        tstops = _tstops_callable
    end

    opts = OrdinaryDiffEqCore.DEOptions(
        maxiters, save_everystep,
        adaptive, abstol_internal,
        reltol_internal,
        QT(gamma),
        QT(qmax), QT(qmin),
        QT(qsteady_max), QT(qsteady_min),
        QT(qoldinit),
        QT(failfactor),
        tType(dtmax), tType(dtmin),
        controller,
        internalnorm, LinearAlgebra.opnorm,
        save_idxs,
        tstops_internal, saveat_internal,
        d_discontinuities_internal,
        tstops, saveat, d_discontinuities,
        userdata,
        progress, progress_steps,
        progress_name, progress_message, progress_id,
        timeseries_errors, dense_errors,
        convert.(uBottomEltypeNoUnits, delta),
        dense, save_on, save_start,
        save_end, save_noise, save_discretes, save_end_user,
        callbacks_internal, isoutofdomain, unstable_check,
        verbose_internal, calck, force_dtmin,
        advance_to_tstop, stop_at_next_tstop
    )

    stats = DiffEqBase.Stats(0)
    if alg isa Union{
            StochasticDiffEqCompositeAlgorithm,
            StochasticDiffEqRODECompositeAlgorithm,
        }
        sol = DiffEqBase.build_solution(
            prob, alg, ts, timeseries, W = W,
            stats = stats, saved_subsystem = saved_subsystem,
            calculate_error = false, alg_choice = alg_choice,
            interp = id, dense = dense, seed = _seed
        )
    else
        sol = DiffEqBase.build_solution(
            prob, alg, ts, timeseries, W = W,
            stats = stats, saved_subsystem = saved_subsystem,
            calculate_error = false,
            interp = id, dense = dense, seed = _seed
        )
    end

    FType = typeof(f)
    CType = typeof(c)
    SolType = typeof(sol)
    cacheType = typeof(cache)

    tprev = t
    dtcache = tType(dt)
    iter = 0
    u_modified = false
    eigen_est = inv(one(tType)) # rate/state = (state/time)/state = 1/t units
    EEst = tTypeNoUnits(1)
    just_hit_tstop = false
    do_error_check = true
    isout = false
    accept_step = false
    force_stepfail = false
    last_stepfail = false
    event_last_time = 0
    vector_event_last_time = 1
    last_event_error = zero(uBottomEltypeNoUnits)
    dtchangeable = true
    q11 = tTypeNoUnits(1)
    success_iter = 0
    q = tTypeNoUnits(1)

    ks = Vector{rateType}(undef, 0)

    integrator = OrdinaryDiffEqCore.ODEIntegrator{
        typeof(alg), isinplace(prob),
        uType, Nothing, tType, typeof(p), typeof(eigen_est), tTypeNoUnits,
        QT, typeof(tdir),
        typeof(ks), typeof(sol), FType, typeof(cache), typeof(opts), Nothing,
        typeof(last_event_error), typeof(callback_cache),
        typeof(initializealg), Nothing, typeof(controller), typeof(_rng),
        typeof(W), typeof(P), tType,
        typeof(noise), CType, typeof(rate_constants),
    }(
        sol, u, nothing, ks, t, tType(dt), f, p,
        uprev, uprev, nothing, tprev,
        alg, dtcache, dtchangeable,
        tType(dt), tdir, eigen_est, EEst,
        QT(qoldinit), q11,
        QT(1), tType(1),
        controller,
        success_iter,
        iter, saveiter, 0, cache,
        callback_cache,
        0, force_stepfail,
        last_stepfail,
        just_hit_tstop, false, zero(t), do_error_check,
        event_last_time,
        vector_event_last_time,
        last_event_error, accept_step,
        isout, false,
        u_modified, true, false,
        opts, stats, initializealg, nothing,
        nothing, nothing, _rng,
        W, P, tType(dt),
        noise, c, rate_constants, q
    )

    if initialize_integrator
        DiffEqBase.initialize_dae!(integrator)
        initialize_callbacks!(integrator, initialize_save)
        initialize!(integrator, integrator.cache)
        save_start &&
            alg isa Union{
            StochasticDiffEqCompositeAlgorithm,
            StochasticDiffEqRODECompositeAlgorithm,
        } &&
            copyat_or_push!(alg_choice, 1, integrator.cache.current)
    end

    # Evaluate callable tstops now that callbacks are initialized and parameters finalized.
    if _tstops_callable !== nothing
        for ts in _tstops_callable(integrator.p, prob.tspan)
            add_tstop!(integrator, ts)
        end
    end

    handle_dt!(integrator)

    ## Modify the first dt for tstops
    modify_dt_for_tstops!(integrator)
    ### Needs to be done before first rand
    integrator.sqdt = integrator.tdir * sqrt(abs(integrator.dt))

    !isnothing(integrator.W) && (integrator.W.dt = integrator.dt)
    !isnothing(integrator.P) && (integrator.P.dt = integrator.dt)

    DiffEqNoiseProcess.setup_next_step!(integrator::SDEIntegrator)

    return integrator
end
end
end

function DiffEqBase.solve!(integrator::SDEIntegrator)
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

    f = integrator.sol.prob.f isa Tuple ? integrator.sol.prob.f[1] : integrator.sol.prob.f

    if DiffEqBase.has_analytic(f)
        DiffEqBase.calculate_solution_errors!(
            integrator.sol; timeseries_errors = integrator.opts.timeseries_errors,
            dense_errors = integrator.opts.dense_errors
        )
    end
    if integrator.sol.retcode != ReturnCode.Default
        return integrator.sol
    end
    return integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, ReturnCode.Success)
end

# Helpers

function handle_dt!(integrator)
    return if iszero(integrator.dt) && integrator.opts.adaptive
        auto_dt_reset!(integrator)
        if sign(integrator.dt) != integrator.tdir && !iszero(integrator.dt) &&
                !isnan(integrator.dt)
            error("Automatic dt setting has the wrong sign. Exiting. Please report this error.")
        end
        if isnan(integrator.dt)
            @SciMLMessage(
                "Automatic dt set the starting dt as NaN, causing instability.",
                integrator.opts.verbose, :dt_NaN
            )
        end
    elseif integrator.opts.adaptive && integrator.dt > zero(integrator.dt) &&
            integrator.tdir < 0
        integrator.dt *= integrator.tdir # Allow positive dt, but auto-convert
    end
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

        if initialize_save &&
                (
                any((c) -> c.save_positions[2], callbacks.discrete_callbacks) ||
                    any((c) -> c.save_positions[2], callbacks.continuous_callbacks)
            )
            savevalues!(integrator, true)
        end
    end

    # reset this as it is now handled so the integrators should proceed as normal
    integrator.u_modified = false

    return if initialize_save
        SciMLBase.save_discretes_if_enabled!(integrator, integrator.opts.callback; skip_duplicates = true)
    end
end

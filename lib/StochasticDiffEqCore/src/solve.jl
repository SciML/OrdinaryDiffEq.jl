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
        if rng isa Random.TaskLocalRNG
            # TaskLocalRNG is the ensemble-layer default, not an explicit
            # per-trajectory RNG. If the user set an explicit seed (via the
            # `seed` kwarg or `prob.seed`, typically from `remake(prob, seed=…)`
            # inside `prob_func`), respect it over the TaskLocalRNG so
            # `remake(prob, seed=s)` inside an ensemble is actually reproducible.
            if !iszero(seed)
                return Random.Xoshiro(seed), seed, false
            elseif prob isa DiffEqBase.AbstractRODEProblem && !iszero(prob.seed)
                return Random.Xoshiro(prob.seed), prob.seed, false
            end
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

"""
    _z_prototype(alg, rand_prototype, iip::Bool) -> rand_prototype2

Compute the Z process prototype for algorithms that need an extra Brownian process.
Default: use rand_prototype itself as Z prototype (same shape as W).
Solver subpackages override this for algorithms with special Z requirements
(e.g., PL1WM, RKMilGeneral, W2Ito1).
"""
function _z_prototype(alg, rand_prototype, iip::Bool)
    return rand_prototype
end

function DiffEqBase.__init(
        _prob::Union{DiffEqBase.AbstractRODEProblem, JumpProblem},
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm};
        kwargs...
    )
    return _sde_init(_prob, alg; kwargs...)
end

function _sde_init(
        _prob::Union{DiffEqBase.AbstractRODEProblem, JumpProblem},
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
        dense = false,
        calck = false,
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
        verbose = Standard(), force_dtmin = false,
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
    # NOTE: JumpProblem kwargs merge is already done by init_call / __solve
    # before reaching _sde_init. DO NOT call merge_problem_kwargs again here.

    is_sde = _prob isa SDEProblem

    # ── Alias resolution (SDE/RODE-specific specifier types) ─────────────
    if alias isa Bool
        aliases = is_sde ? SciMLBase.SDEAliasSpecifier(; alias) :
            SciMLBase.RODEAliasSpecifier(; alias)
    elseif alias isa SciMLBase.SDEAliasSpecifier
        aliases = alias
    elseif alias isa SciMLBase.RODEAliasSpecifier
        aliases = alias
    else
        aliases = is_sde ? SciMLBase.SDEAliasSpecifier() :
            SciMLBase.RODEAliasSpecifier()
    end

    prob = concrete_prob(_prob)

    # ── RNG resolution ───────────────────────────────────────────────────
    _rng, _seed, _rng_provided = _resolve_rng(rng, seed, prob)

    # ── JumpProblem reset ────────────────────────────────────────────────
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
    prob = concrete_prob(_prob)

    # ── SDE-specific validation ──────────────────────────────────────────
    if typeof(prob.f) <: Tuple
        if any(mm != I for mm in prob.f.mass_matrix)
            error("This solver is not able to use mass matrices.")
        end
    elseif prob isa DiffEqBase.AbstractRODEProblem && prob.f.mass_matrix != I &&
            !alg_mass_matrix_compatible(alg)
        error("This solver is not able to use mass matrices.")
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

    # ── StochasticCompositeAlgorithm AutoSwitch → AutoSwitchCache ────────
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

    # ── f/p/u aliasing (needed before cache construction) ────────────────
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

    if prob.u0 isa Tuple
        u = ArrayPartition(prob.u0, Val{true})
    end

    # ── Type computations (needed for cache) ─────────────────────────────
    tType = eltype(prob.tspan)
    uType = typeof(u)
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)
    uEltypeNoUnits = recursive_unitless_eltype(u)
    tTypeNoUnits = typeof(one(tType))
    noise = prob isa DiffEqBase.AbstractRODEProblem ? prob.noise : nothing
    tspan = prob.tspan
    t = tspan[1]

    # ── SDE-specific abstol/reltol defaults (1//10^2 vs ODE's 1//10^6) ──
    if abstol === nothing
        if uBottomEltypeNoUnits === uBottomEltype && !(uBottomEltype <: Integer)
            abstol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        elseif uBottomEltype <: Integer
            abstol = real(oneunit(uBottomEltype) * 1 // 10^2)
        else
            abstol = 0
        end
    end

    if reltol === nothing
        if uBottomEltypeNoUnits === uBottomEltype && !(uBottomEltype <: Integer)
            reltol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1 // 10^2))
        elseif uBottomEltype <: Integer
            reltol = real(oneunit(uBottomEltype) * 1 // 10^2)
        else
            reltol = 0
        end
    end

    # ── rate_prototype / noise_rate_prototype (needed for cache) ─────────
    if isinplace(prob) && u isa AbstractArray && eltype(u) <: Number &&
            uBottomEltypeNoUnits == uBottomEltype
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
        else
            rate_prototype = u / oneunit(tType)
        end
    end
    rateType = typeof(rate_prototype)

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

    # ── uprev (needed for cache) ─────────────────────────────────────────
    uprev = recursivecopy(u)

    # ── rand_prototype (needed for noise creation) ───────────────────────
    if !(uType <: AbstractArray)
        rand_prototype = zero(u ./ u)
    else
        randElType = uBottomEltypeNoUnits
        if prob.f isa DynamicalSDEFunction
            rand_prototype = copy(noise_rate_prototype)
        elseif is_diagonal_noise(prob)
            if u isa SArray
                rand_prototype = zero(u)
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
    end

    # ── Callback merging for JumpProblem ─────────────────────────────────
    _callback = _prob isa JumpProblem ?
        CallbackSet(callback, _prob.jump_callback) : callback

    # ── Noise creation (WienerProcess / user noise handling) ─────────────
    if prob isa DiffEqBase.AbstractRODEProblem && prob.noise === nothing
        rswm = isadaptive(alg) ? RSWM(adaptivealg = :RSwM3) : RSWM(adaptivealg = :RSwM1)
        if isinplace(prob)
            if alg_needs_extra_process(alg)
                rand_prototype2 = _z_prototype(alg, rand_prototype, true)
                W = WienerProcess!(
                    t, rand_prototype, rand_prototype2,
                    save_everystep = save_noise,
                    rng = _rng
                )
            else
                W = WienerProcess!(
                    t, rand_prototype,
                    save_everystep = save_noise,
                    rng = _rng
                )
            end
        else
            if alg_needs_extra_process(alg)
                rand_prototype2 = _z_prototype(alg, rand_prototype, false)
                W = WienerProcess(
                    t, rand_prototype, rand_prototype2,
                    save_everystep = save_noise,
                    rng = _rng
                )
            else
                W = WienerProcess(
                    t, rand_prototype,
                    save_everystep = save_noise,
                    rng = _rng
                )
            end
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
            if !_rng_provided && W isa Union{NoiseProcess, NoiseTransport} && W.reseed
                Random.seed!(W.rng, _seed)
            end
            if W.curt != t
                reinit!(W, t, t0 = t)
            end
        elseif W.curt != t
            error("Starting time in the noise process is not the starting time of the simulation. The noise process should be re-initialized for repeated use")
        end
    else
        @assert _prob isa JumpProblem
        W = nothing
    end

    # ── CompoundPoissonProcess for tau-leaping ───────────────────────────
    if _prob isa JumpProblem && _prob.regular_jump !== nothing
        if !isnothing(_prob.regular_jump.mark_dist) == nothing
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

    # ── dW/dZ extraction, verbose conversion, alg_cache ──────────────────
    dW, dZ = isnothing(W) ? (nothing, nothing) : (W.dW, W.dZ)

    verbose_internal = _process_verbose_param(verbose)

    cache = alg_cache(
        alg, prob, u, dW, dZ, p, rate_prototype, noise_rate_prototype,
        jump_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, f, t, dt, Val{isinplace(_prob)}, verbose_internal
    )

    # ── Delegate to ODE's _ode_init directly ─────────────────────────────
    ode_alias = ODEAliasSpecifier(alias_u0 = true, alias_f = true, alias_p = true)

    tType = eltype(prob.tspan)

    integrator = OrdinaryDiffEqCore._ode_init(
        prob, alg;
        # Pre-built objects
        _cache = cache,
        _u = u,
        _uprev = uprev,
        W = W, P = P,
        sqdt = tType(dt),
        noise = noise, c = c, rate_constants = rate_constants,
        seed = _seed,
        rng = _rng,
        controller,
        # SDE-specific defaults (passed through from SDE's kwargs)
        saveat, tstops, d_discontinuities,
        save_idxs, save_everystep, save_noise, save_on,
        save_start, save_end, callback = _callback,
        dense, calck, dt = iszero(dt) ? nothing : dt, adaptive,
        abstol, reltol,
        fullnormalize, failfactor,
        delta = convert.(uBottomEltypeNoUnits, delta),
        maxiters, dtmax, dtmin,
        internalnorm, isoutofdomain, unstable_check,
        verbose, force_dtmin,
        timeseries_errors, dense_errors,
        advance_to_tstop, stop_at_next_tstop,
        initialize_save, progress, progress_steps,
        progress_name, progress_message, progress_id,
        userdata, save_discretes, initialize_integrator,
        alias = ode_alias,
        initializealg,
        kwargs...
    )

    # ── Post-delegation: noise process next-step setup ───────────────────
    DiffEqNoiseProcess.setup_next_step!(integrator::SDEIntegrator)

    return integrator
end

# solve! is now provided by OrdinaryDiffEqCore (SciMLBase.solve!(::ODEIntegrator))
# handle_dt! and initialize_callbacks! are now provided by OrdinaryDiffEqCore

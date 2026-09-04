using Statistics

# Errors recorded for a tolerance whose solve failed or was cut off by `timeout`.
_nan_errors() = Dict{Symbol, Float64}(:l∞ => NaN, :L2 => NaN, :final => NaN, :l2 => NaN, :L∞ => NaN)

_setup_tags(setup) = Vector{Symbol}(get(setup, :tags, Symbol[]))
_setup_name(setup) = get(setup, :name, _default_name(setup[:alg]))

_preset_algorithm_tags(alg) = Symbol[]

const _FAMILY_TAGS = Set(
    [
        :adams, :bdf, :esdirk, :exponential, :explicit_rk, :extrapolation, :firk, :imex,
        :low_storage_rk, :multistep, :partitioned, :rk, :rkn, :rosenbrock, :sdirk,
        :ssprk, :stabilized, :stabilized_rk, :symplectic,
    ]
)
const _TRAIT_TAGS = Set(
    [
        :adaptive, :composite, :explicit, :fixed_step, :implicit, :split, :variable_order,
    ]
)
const _ROLE_TAGS = Set([:recommended, :reference])
const _PROVIDER_TAGS = Set([:external, :lsoda, :odeinterface, :ordinarydiffeq, :sundials])
const _VARIANT_TAGS = Set(
    [
        :dae_residual, :dense, :dense_output, :dynamic, :endpoint, :mass_matrix,
        :matrix_free, :sparse, :static, :timeseries,
    ]
)
const _DOMAIN_TAGS = Set([:bv, :dae, :dde, :nonstiff, :ode, :pde, :rode, :sde, :stiff])

"""
    tag_kind(tag::Symbol) -> Symbol

Return the semantic kind of a benchmark tag: `:family`, `:trait`, `:role`, `:provider`,
`:variant`, `:domain`, or `:unknown`. Order tags such as `:order_5` are traits, and
automatic-differentiation tags such as `:autodiff_forwarddiff` are variants.

[`autoplot`](@ref) uses only `:family` tags for automatic family discovery. Pass its
`families` keyword for one-off custom families. Packages defining a stable family tag
can make it discoverable with a specialization such as
`DiffEqDevTools.tag_kind(::Val{:my_family}) = :family`.
"""
tag_kind(tag::Symbol) = tag_kind(Val(tag))
function tag_kind(::Val{tag}) where {tag}
    tag in _FAMILY_TAGS && return :family
    (tag in _TRAIT_TAGS || startswith(string(tag), "order_")) && return :trait
    tag in _ROLE_TAGS && return :role
    tag in _PROVIDER_TAGS && return :provider
    (tag in _VARIANT_TAGS || startswith(string(tag), "autodiff_")) && return :variant
    tag in _DOMAIN_TAGS && return :domain
    return :unknown
end

function _known_alg_order(alg)
    applicable(SciMLBase.alg_order, alg) || return nothing
    try
        return SciMLBase.alg_order(alg)
    catch err
        if err isa ErrorException && err.msg == "Order is not defined for this algorithm"
            return nothing
        end
        rethrow()
    end
end

_order_tag(order) = Symbol("order_", replace(string(order), "//" => "_", "." => "_"))

"""
    auto_tags(alg) -> Vector{Symbol}

Return preset metadata tags derived from the public traits and type hierarchy of `alg`.
Differential-equation algorithms receive an order tag such as `:order_5` (or
`:order_3_2` for order `3 // 2`) when their order is defined, plus `:adaptive` or
`:fixed_step`. When OrdinaryDiffEqCore is loaded, the result also describes known
implicitness, families, and structural traits such as `:rosenbrock`, `:firk`, `:sdirk`,
`:split`, and `:multistep`.

Algorithms without applicable traits return an empty vector. [`WorkPrecision`](@ref)
and [`WorkPrecisionSet`](@ref) merge these presets with explicit tags by default.
"""
auto_tags(alg) = Symbol[]
function auto_tags(alg::AbstractDEAlgorithm)
    presets = _preset_algorithm_tags(alg)
    tags = Symbol[]
    order = :variable_order in presets ? nothing : _known_alg_order(alg)
    order === nothing || push!(tags, _order_tag(order))
    push!(tags, SciMLBase.isadaptive(alg) ? :adaptive : :fixed_step)
    append!(tags, presets)
    return unique!(tags)
end

function _combined_tags(alg, tags, include_auto_tags)
    combined = include_auto_tags ? auto_tags(alg) : Symbol[]
    append!(combined, tags)
    return unique!(combined)
end

function _timed_out(timeout, solve_time, name, abstol)
    (timeout === nothing || solve_time <= timeout) && return false
    @warn "$name exceeded the $(timeout)s timeout at abstol=$abstol " *
        "($(round(solve_time, sigdigits = 3))s); recording this point as NaN."
    return true
end

# Build the `errors` StructArray from the per-tolerance error `Dict`s.
# A failed tolerance records a different key set than a successful one and `Dict`
# iteration order is unspecified, so `NamedTuple.(dicts)` can yield differently-typed
# NamedTuples whose element type widens to bare `NamedTuple` — which StructArrays 0.7
# rejects. Use the union of keys in a fixed order and fill the gaps with NaN.
function _dicts_to_structarray(dicts)
    ks = Tuple(sort!(unique!(mapreduce(d -> collect(keys(d)), vcat, dicts))))
    V = mapreduce(valtype, promote_type, dicts)
    T = all(d -> length(d) == length(ks), dicts) ? V : promote_type(V, Float64)
    NT = NamedTuple{ks, NTuple{length(ks), T}}
    return StructArray([NT(Tuple(get(d, k, NaN) for k in ks)) for d in dicts])
end

# Default names of algorithms:
# Workaround for `MethodOfSteps` algorithms, otherwise they are all called "MethodOfSteps"
# Ideally this would be a trait (in SciMLBase?), so packages could implement it
_default_name(alg) = __default_name(alg)
function _default_name(alg::AbstractDDEAlgorithm)
    return isdefined(alg, :alg) ? _default_name(alg.alg) : __default_name(alg)
end
__default_name(alg) = string(nameof(typeof(alg)))

## Shootouts

"""
    Shootout(
        prob, setups; appxsol = nothing, names = nothing,
        error_estimate = :final, numruns = 20, seconds = 2, kwargs...
    )

Benchmark multiple solver configurations on one problem. Each entry of `setups` is a
dictionary containing an `:alg` and any solver-specific keyword arguments. The result
stores the solutions, errors, timings, efficiencies, and the configuration with the
highest efficiency, where efficiency is `1 / (error * time)`.

Use `appxsol` as a numerical reference when the problem has no analytic solution.
Additional keyword arguments are forwarded to every solve.
"""
mutable struct Shootout
    setups::Vector{Dict{Symbol, Any}}
    times::Any #::Vector{Float64}
    errors::Any #::Vector{uType}
    effs::Any #::Vector{Float64} # Efficiencies
    effratios::Any #::Matrix{uEltype}
    solutions::Any
    names::Vector{String}
    N::Int
    bestidx::Int
    winner::String
end

"""
    ShootoutSet(
        probs, setups; probaux = nothing, names = nothing,
        print_names = false, kwargs...
    )

Run a [`Shootout`](@ref) for every problem in `probs`. `probaux` may contain one
dictionary of per-problem keyword arguments for each problem; other keyword arguments
are shared by every shootout.
"""
mutable struct ShootoutSet
    shootouts::Vector{Shootout}
    probs::Any #::Vector{AbstractDEProblem}
    probaux::Any #::Vector{Dict{Symbol,Any}}
    N::Int
    winners::Vector{String}
end

function ode_shootout(args...; kwargs...)
    @warn("ode_shootout is deprecated. Use ShootOut instead")
    return ShootOut(args...; kwargs...)
end

function Shootout(
        prob, setups; appxsol = nothing, names = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, kwargs...
    )
    N = length(setups)
    @assert names === nothing || length(setups) == length(names)
    errors = Vector{Float64}(undef, N)
    solutions = Vector{Any}(undef, N)
    effs = Vector{Float64}(undef, N)
    times = Vector{Float64}(undef, N)
    effratios = Matrix{Float64}(undef, N, N)
    timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
    dense_errors = error_estimate ∈ DENSE_ERRORS
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
    end
    for i in eachindex(setups)
        sol = solve(
            prob, setups[i][:alg]; timeseries_errors,
            dense_errors, kwargs..., setups[i]...
        )

        if :prob_choice ∈ keys(setups[i])
            cur_appxsol = appxsol[setups[i][:prob_choice]]
        elseif prob isa AbstractArray
            cur_appxsol = appxsol[1]
        else
            cur_appxsol = appxsol
        end

        if cur_appxsol !== nothing
            errsol = appxtrue(sol, cur_appxsol)
            errors[i] = errsol.errors[error_estimate]
            solutions[i] = errsol
        else
            errors[i] = sol.errors[error_estimate]
            solutions[i] = sol
        end

        if haskey(setups[i], :prob_choice)
            _prob = prob[setups[i][:prob_choice]]
        elseif prob isa AbstractArray
            _prob = prob[1]
        else
            _prob = prob
        end

        benchmark_f = let _prob = _prob, alg = setups[i][:alg], sol = sol, kwargs = kwargs
            () -> @elapsed solve(
                _prob, alg, sol.u, sol.t, sol.k;
                timeseries_errors = false,
                dense_errors = false, kwargs...
            )
        end
        benchmark_f() # pre-compile

        b_t = benchmark_f()
        if b_t > seconds
            times[i] = b_t
        else
            times[i] = mapreduce(i -> benchmark_f(), min, 2:numruns; init = b_t)
        end

        effs[i] = 1 / (errors[i] * times[i])
    end
    for j in 1:N, i in 1:N
        effratios[i, j] = effs[i] / effs[j]
    end
    bestidx = findall((y) -> y == maximum(effs), effs)[1]
    winner = names[bestidx]
    return Shootout(
        setups, times, errors, effs, effratios, solutions, names, N, bestidx,
        winner
    )
end

function ode_shootoutset(args...; kwargs...)
    @warn("ode_shootoutset is deprecated. Use ShootoutSet instead")
    return ShootoutSet(args...; kwargs...)
end

function ShootoutSet(
        probs, setups; probaux = nothing,
        names = nothing, print_names = false, kwargs...
    )
    N = length(probs)
    @assert names === nothing || length(setups) == length(names)
    shootouts = Vector{Shootout}(undef, N)
    winners = Vector{String}(undef, N)
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
    end
    if probaux === nothing
        probaux = Vector{Dict{Symbol, Any}}(undef, N)
        for i in 1:N
            probaux[i] = Dict{Symbol, Any}()
        end
    end
    for i in eachindex(probs)
        print_names && println(names[i])
        shootouts[i] = Shootout(probs[i], setups; names, kwargs..., probaux[i]...)
        winners[i] = shootouts[i].winner
    end
    return ShootoutSet(shootouts, probs, probaux, N, winners)
end

Base.length(shoot::Shootout) = shoot.N
Base.size(shoot::Shootout) = length(shoot)
Base.getindex(shoot::Shootout, i::Int) = shoot.effs[i]
Base.getindex(shoot::Shootout, ::Colon) = shoot.effs
Base.firstindex(shoot::Shootout) = 1
Base.lastindex(shoot::Shootout) = lastindex(shoot.effs)

function Base.show(io::IO, shoot::Shootout)
    println(io, "Winner: $(shoot.winner)")
    return println(io, "EffRatios: $(shoot.effratios[shoot.bestidx, :])")
end

Base.length(set::ShootoutSet) = set.N
Base.size(set::ShootoutSet) = length(set)
Base.getindex(set::ShootoutSet, i::Int) = set.shootouts[i]
Base.getindex(set::ShootoutSet, ::Colon) = set.shootouts
Base.show(io::IO, set::ShootoutSet) = print(io, "ShootoutSet of $(set.N) shootouts ")
Base.firstindex(shoot::ShootoutSet) = 1
Base.lastindex(shoot::ShootoutSet) = lastindex(shoot.shootouts)

## WorkPrecisions

"""
    WorkPrecision(
        prob, alg, abstols, reltols, dts = nothing;
        name = nothing, appxsol = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, tags = Symbol[], auto_tags = true,
        timeout = nothing, kwargs...
    )

Measure error and execution time for `alg` at corresponding absolute and relative
tolerances. When `dts` is provided, its entries select a fixed time step for each
tolerance pair. The result stores the measured errors, timings, solver statistics, and
inputs for plotting a work-precision diagram.

Use `appxsol` as a numerical reference when the problem has no analytic solution.
`tags` attaches metadata symbols used by [`filter_by_tags`](@ref) and the plot recipe.
They are merged with [`auto_tags(alg)`](@ref) unless `auto_tags = false`.
`timeout` gives a per-tolerance wall-clock budget in seconds; see
[`WorkPrecisionSet`](@ref). Additional keyword arguments are forwarded to `solve`.
"""
mutable struct WorkPrecision
    prob::Any
    abstols::Any
    reltols::Any
    errors::Any
    times::Any
    dts::Any
    stats::Any
    name::Any
    error_estimate::Any
    N::Int
    tags::Vector{Symbol}
end

function WorkPrecision(
        prob, abstols, reltols, errors, times, dts, stats, name, error_estimate, N
    )
    return WorkPrecision(
        prob, abstols, reltols, errors, times, dts, stats, name, error_estimate, N,
        Symbol[]
    )
end

"""
    WorkPrecisionSet(prob, abstols, reltols, setups; error_estimates = nothing,
        timeout = nothing, kwargs...)

Build one [`WorkPrecision`](@ref) result for each solver configuration in `setups` so
their work-precision curves can be compared. Each setup is a dictionary containing an
`:alg` and may override the shared tolerances or fixed step sizes with `:abstols`,
`:reltols`, or `:dts`, or its legend entry with `:name`. Preset tags from
[`auto_tags`](@ref) are attached by default and merged with any `:tags` entry. Set
`:auto_tags => false` in a setup to use only its explicit tags. The resulting metadata
is used by [`filter_by_tags`](@ref), [`best_of_families`](@ref), [`autoplot`](@ref),
and the plot recipe to build family and cross-family comparisons.

`error_estimates` requests several error metrics from a single run (for example
`[:final, :l2, :L2]`), so plots for each metric can be drawn without re-solving; the
computed metrics are reported by [`available_errors`](@ref). `error_estimate` remains
the metric the plot recipe defaults to.

`timeout` is a per-tolerance wall-clock budget in seconds. A solve is never
interrupted, but once one exceeds the budget its repeated timing runs are skipped and
its error and time are recorded as `NaN`, which the plot recipe drops.
"""
mutable struct WorkPrecisionSet
    wps::Vector{WorkPrecision}
    N::Int
    abstols::Any
    reltols::Any
    prob::Any
    setups::Any
    names::Any
    error_estimate::Any
    numruns::Any
    active_error_estimates::Vector{Symbol}
end

function WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate, numruns
    )
    return WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate, numruns,
        Symbol[error_estimate]
    )
end

"""
    available_errors(wp_set::WorkPrecisionSet) -> Vector{Symbol}

Return the error estimates computed for `wp_set`, i.e. the `error_estimates` requested
from [`WorkPrecisionSet`](@ref) or, when none were, the single `error_estimate`. Each
one is a valid `x` for the plot recipe.
"""
available_errors(wp_set::WorkPrecisionSet) = wp_set.active_error_estimates

_needs_timeseries(error_estimates) = any(e -> e ∈ TIMESERIES_ERRORS, error_estimates)
_needs_dense(error_estimates) = any(e -> e ∈ DENSE_ERRORS, error_estimates)

function WorkPrecision(
        prob, alg, abstols, reltols, dts = nothing;
        name = nothing, appxsol = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, timeout = nothing,
        timeseries_errors::Union{Bool, Nothing} = nothing,
        dense_errors::Union{Bool, Nothing} = nothing,
        tags::Vector{Symbol} = Symbol[], auto_tags::Bool = true, kwargs...
    )
    N = length(abstols)
    errors = Vector{Dict{Symbol, Float64}}(undef, N)
    times = Vector{Float64}(undef, N)
    stats = Vector{Any}(undef, N)
    if name === nothing
        name = "WP-Alg"
    end

    if haskey(kwargs, :prob_choice)
        _prob = prob[kwargs[:prob_choice]]
    elseif prob isa AbstractArray
        _prob = prob[1]
    else
        _prob = prob
    end

    let _prob = _prob
        timeseries_errors = something(timeseries_errors, error_estimate ∈ TIMESERIES_ERRORS)
        dense_errors = something(dense_errors, error_estimate ∈ DENSE_ERRORS)
        for i in 1:N
            t_start = time()
            if dts === nothing
                sol = solve(
                    _prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], timeseries_errors,
                    dense_errors
                )
            else
                sol = solve(
                    _prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], dt = dts[i],
                    timeseries_errors,
                    dense_errors
                )
            end

            stats[i] = sol.stats

            if _timed_out(timeout, time() - t_start, name, abstols[i])
                errors[i] = _nan_errors()
                times[i] = NaN
                continue
            end

            if SciMLBase.successful_retcode(sol)
                if haskey(kwargs, :prob_choice)
                    cur_appxsol = appxsol[kwargs[:prob_choice]]
                elseif prob isa AbstractArray
                    cur_appxsol = appxsol[1]
                else
                    cur_appxsol = appxsol
                end

                if cur_appxsol !== nothing
                    errsol = appxtrue(sol, cur_appxsol)
                    errors[i] = Dict{Symbol, Float64}()
                    for err in keys(errsol.errors)
                        errors[i][err] = mean(errsol.errors[err])
                    end
                else
                    errors[i] = Dict{Symbol, Float64}()
                    for err in keys(sol.errors)
                        errors[i][err] = mean(sol.errors[err])
                    end
                end

                benchmark_f = let dts = dts, _prob = _prob, alg = alg, sol = sol,
                        abstols = abstols, reltols = reltols, kwargs = kwargs

                    if dts === nothing
                        if _prob isa DAEProblem
                            () -> @elapsed solve(
                                _prob, alg, sol.u, sol.t;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        else
                            () -> @elapsed solve(
                                _prob, alg, sol.u, sol.t, sol.k;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        end
                    else
                        if _prob isa DAEProblem
                            () -> @elapsed solve(
                                _prob, alg, sol.u, sol.t;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        else
                            () -> @elapsed solve(
                                _prob, alg, sol.u, sol.t, sol.k;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        end
                    end
                end
                benchmark_f() # pre-compile

                b_t = benchmark_f()
                if b_t > seconds
                    times[i] = b_t
                else
                    times[i] = mapreduce(i -> benchmark_f(), min, 2:numruns; init = b_t)
                end
            else
                # Unsuccessful retcode, give NaN time
                errors[i] = _nan_errors()
                times[i] = NaN
            end
        end
    end
    return WorkPrecision(
        prob, abstols, reltols, _dicts_to_structarray(errors),
        times, dts, stats, name, error_estimate, N,
        _combined_tags(alg, tags, auto_tags)
    )
end

# Work precision information for a BVP
function WorkPrecision(
        prob::AbstractBVProblem, alg, abstols, reltols, dts = nothing;
        name = nothing, appxsol = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, timeout = nothing,
        timeseries_errors::Union{Bool, Nothing} = nothing,
        dense_errors::Union{Bool, Nothing} = nothing,
        tags::Vector{Symbol} = Symbol[], auto_tags::Bool = true, kwargs...
    )
    N = length(abstols)
    errors = Vector{Dict{Symbol, Float64}}(undef, N)
    times = Vector{Float64}(undef, N)
    stats = Vector{Any}(undef, N)
    if name === nothing
        name = "WP-Alg"
    end

    if haskey(kwargs, :prob_choice)
        _prob = prob[kwargs[:prob_choice]]
    elseif prob isa AbstractArray
        _prob = prob[1]
    else
        _prob = prob
    end

    let _prob = _prob
        timeseries_errors = something(timeseries_errors, error_estimate ∈ TIMESERIES_ERRORS)
        dense_errors = something(dense_errors, error_estimate ∈ DENSE_ERRORS)
        for i in 1:N
            t_start = time()
            if dts === nothing
                sol = solve(
                    _prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], timeseries_errors,
                    dense_errors
                )
            else
                sol = solve(
                    _prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], dt = dts[i],
                    timeseries_errors,
                    dense_errors
                )
            end

            stats[i] = sol.stats

            if _timed_out(timeout, time() - t_start, name, abstols[i])
                errors[i] = _nan_errors()
                times[i] = NaN
                continue
            end

            if SciMLBase.successful_retcode(sol)
                if haskey(kwargs, :prob_choice)
                    cur_appxsol = appxsol[kwargs[:prob_choice]]
                elseif prob isa AbstractArray
                    cur_appxsol = appxsol[1]
                else
                    cur_appxsol = appxsol
                end

                if cur_appxsol !== nothing
                    errsol = appxtrue(sol, cur_appxsol)
                    errors[i] = Dict{Symbol, Float64}()
                    for err in keys(errsol.errors)
                        errors[i][err] = mean(errsol.errors[err])
                    end
                else
                    errors[i] = Dict{Symbol, Float64}()
                    for err in keys(errsol.errors)
                        errors[i][err] = mean(errsol.errors[err])
                    end
                end

                benchmark_f = let dts = dts, _prob = _prob, alg = alg, sol = sol,
                        abstols = abstols, reltols = reltols, kwargs = kwargs

                    if dts === nothing
                        if _prob isa DAEProblem
                            () -> @elapsed solve(
                                _prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        else
                            () -> @elapsed solve(
                                _prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        end
                    else
                        if _prob isa DAEProblem
                            () -> @elapsed solve(
                                _prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        else
                            () -> @elapsed solve(
                                _prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...
                            )
                        end
                    end
                end
                benchmark_f() # pre-compile

                b_t = benchmark_f()
                if b_t > seconds
                    times[i] = b_t
                else
                    times[i] = mapreduce(i -> benchmark_f(), min, 2:numruns; init = b_t)
                end
            else
                # Unsuccessful retcode, give NaN error and time
                errors[i] = _nan_errors()
                times[i] = NaN
            end
        end
    end
    return WorkPrecision(
        prob, abstols, reltols, _dicts_to_structarray(errors),
        times, dts, stats, name, error_estimate, N,
        _combined_tags(alg, tags, auto_tags)
    )
end

# Work precision information for a nonlinear problem.
function WorkPrecision(
        prob::NonlinearProblem, alg, abstols, reltols, dts = nothing; name = nothing,
        appxsol = nothing, error_estimate = :l2, numruns = 20, seconds = 2,
        timeout = nothing, tags::Vector{Symbol} = Symbol[], auto_tags::Bool = true,
        kwargs...
    )
    N = length(abstols)
    errors = Vector{Dict{Symbol, Float64}}(undef, N)
    times = Vector{Float64}(undef, N)
    stats = Vector{Any}(undef, N)
    if name === nothing
        name = "WP-Alg"
    end

    if haskey(kwargs, :prob_choice)
        _prob = prob[kwargs[:prob_choice]]
    elseif prob isa AbstractArray
        _prob = prob[1]
    else
        _prob = prob
    end

    let _prob = _prob
        for i in 1:N
            t_start = time()
            sol = solve(_prob, alg; kwargs..., abstol = abstols[i], reltol = reltols[i])

            stats[i] = sol.stats

            if _timed_out(timeout, time() - t_start, name, abstols[i])
                errors[i] = Dict{Symbol, Float64}(error_estimate => NaN)
                times[i] = NaN
                continue
            end

            err = appxsol === nothing ? sol.resid : (sol .- appxsol)
            if error_estimate == :l2
                errors[i] = Dict(error_estimate => norm(err, 2))
            elseif error_estimate == :l∞ || error_estimate == :linf
                errors[i] = Dict(error_estimate => norm(err, Inf))
            else
                error("Unsupported norm used: $(error_estimate).")
            end

            benchmark_f = let dts = dts, _prob = _prob, alg = alg, sol = sol,
                    abstols = abstols, reltols = reltols, kwargs = kwargs

                () -> @elapsed solve(
                    _prob, alg;
                    abstol = abstols[i],
                    reltol = reltols[i],
                    kwargs...
                )
            end
            benchmark_f() # pre-compile

            b_t = benchmark_f()
            if b_t > seconds
                times[i] = b_t
            else
                times[i] = mapreduce(i -> benchmark_f(), min, 2:numruns; init = b_t)
            end
        end
    end

    return WorkPrecision(
        prob, abstols, reltols, _dicts_to_structarray(errors),
        times, dts, stats, name, error_estimate, N,
        _combined_tags(alg, tags, auto_tags)
    )
end

function WorkPrecisionSet(
        prob,
        abstols, reltols, setups;
        print_names = false, names = nothing, appxsol = nothing,
        error_estimate = :final, error_estimates = nothing,
        test_dt = nothing, timeout = nothing, kwargs...
    )
    N = length(setups)
    @assert names === nothing || length(setups) == length(names)
    wps = Vector{WorkPrecision}(undef, N)
    if names === nothing
        names = [_setup_name(setup) for setup in setups]
    end

    active = error_estimates === nothing ? Symbol[error_estimate] : collect(error_estimates)
    timeseries_errors = error_estimates === nothing ? nothing : _needs_timeseries(active)
    dense_errors = error_estimates === nothing ? nothing : _needs_dense(active)

    for i in 1:N
        print_names && println(names[i])
        _abstols = get(setups[i], :abstols, abstols)
        _reltols = get(setups[i], :reltols, reltols)
        _dts = get(setups[i], :dts, nothing)
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[i])

        wps[i] = WorkPrecision(
            prob, setups[i][:alg], _abstols, _reltols, _dts;
            appxsol,
            error_estimate,
            timeout,
            timeseries_errors,
            dense_errors,
            tags = _setup_tags(setups[i]),
            auto_tags = get(setups[i], :auto_tags, true),
            name = names[i], kwargs..., filtered_setup...
        )
    end
    return WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate,
        nothing, active
    )
end

@def error_calculation begin
    if !SciMLBase.has_analytic(prob.f)
        t = prob.tspan[1]:test_dt:prob.tspan[2]
        brownian_values = cumsum(
            [
                [zeros(size(prob.u0))];
                [
                    sqrt(test_dt) * randn(size(prob.u0))
                        for i in 1:(length(t) - 1)
                ]
            ]
        )
        brownian_values2 = cumsum(
            [
                [zeros(size(prob.u0))];
                [
                    sqrt(test_dt) * randn(size(prob.u0))
                        for i in 1:(length(t) - 1)
                ]
            ]
        )
        np = NoiseGrid(t, brownian_values, brownian_values2)
        _prob = remake(prob, noise = np)
        true_sol = solve(_prob, appxsol_setup[:alg]; kwargs..., appxsol_setup...)
    else
        _prob = prob
    end

    # Get a cache
    _abstols = get(setups[1], :abstols, abstols)
    _reltols = get(setups[1], :reltols, reltols)
    _dts = get(setups[1], :dts, zeros(length(_abstols)))
    filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[1])

    sol = solve(
        _prob, setups[1][:alg];
        kwargs..., filtered_setup..., abstol = _abstols[1],
        reltol = _reltols[1], dt = _dts[1],
        timeseries_errors = false,
        dense_errors = false
    )

    for j in 1:M, k in 1:N
        _abstols = get(setups[k], :abstols, abstols)
        _reltols = get(setups[k], :reltols, reltols)
        _dts = get(setups[k], :dts, zeros(length(_abstols)))
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[k])

        sol = solve(
            _prob, setups[k][:alg];
            kwargs..., filtered_setup..., abstol = _abstols[j],
            reltol = _reltols[j], dt = _dts[j],
            timeseries_errors,
            dense_errors
        )
        SciMLBase.has_analytic(prob.f) ? err_sol = sol : err_sol = appxtrue(sol, true_sol)
        tmp_solutions[i, j, k] = err_sol
    end
end

function WorkPrecisionSet(
        prob::AbstractRODEProblem, abstols, reltols, setups,
        test_dt = nothing;
        numruns = 20, numruns_error = 20,
        print_names = false, names = nothing, appxsol_setup = nothing,
        error_estimate = :final, error_estimates = nothing, parallel_type = :none,
        timeout = nothing,
        kwargs...
    )
    @assert names === nothing || length(setups) == length(names)
    timeseries_errors = SciMLBase.has_analytic(prob.f) &&
        error_estimate ∈ TIMESERIES_ERRORS
    weak_timeseries_errors = error_estimate ∈ WEAK_TIMESERIES_ERRORS
    weak_dense_errors = error_estimate ∈ WEAK_DENSE_ERRORS
    dense_errors = SciMLBase.has_analytic(prob.f) && error_estimate ∈ DENSE_ERRORS
    N = length(setups)
    M = length(abstols)
    times = Array{Float64}(undef, M, N)
    tmp_solutions = Array{Any}(undef, numruns_error, M, N)
    if names === nothing
        names = [_setup_name(setup) for setup in setups]
    end
    time_tmp = Vector{Float64}(undef, numruns)

    # First calculate all of the errors
    if parallel_type == :threads
        Threads.@threads for i in 1:numruns_error
            @error_calculation
        end
    elseif parallel_type == :none
        for i in 1:numruns_error
            @info "Error calculation: $i/$numruns_error"
            @error_calculation
        end
    end

    _solutions_k = [
        [EnsembleSolution(tmp_solutions[:, j, k], 0.0, true) for j in 1:M]
            for k in 1:N
    ]
    solutions = [
        [
            SciMLBase.calculate_ensemble_errors(
                sim;
                weak_timeseries_errors,
                weak_dense_errors
            )
                for sim in sol_k
        ] for sol_k in _solutions_k
    ]
    if error_estimate ∈ WEAK_ERRORS
        errors = [[solutions[j][i].weak_errors for i in 1:M] for j in 1:N]
    else
        errors = [[solutions[j][i].error_means for i in 1:M] for j in 1:N]
    end

    local _sol

    # Now time it
    _abstols = [get(setups[k], :abstols, abstols) for k in 1:N]
    _reltols = [get(setups[k], :reltols, reltols) for k in 1:N]
    _dts = [get(setups[k], :dts, zeros(length(_abstols[k]))) for k in 1:N]
    for k in 1:N
        # precompile
        GC.gc()
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[k])

        _sol = solve(
            prob, setups[k][:alg];
            kwargs..., filtered_setup..., abstol = _abstols[k][1],
            reltol = _reltols[k][1], dt = _dts[k][1],
            timeseries_errors = false,
            dense_errors = false
        )
        x = isempty(_sol.t) ? 0 : round(Int, mean(_sol.t) - sum(_sol.t) / length(_sol.t))
        GC.gc()
        for j in 1:M
            timed_out = false
            for i in 1:numruns
                time_tmp[i] = @elapsed sol = solve(
                    prob, setups[k][:alg];
                    kwargs..., filtered_setup...,
                    abstol = _abstols[k][j],
                    reltol = _reltols[k][j], dt = _dts[k][j],
                    timeseries_errors = false,
                    dense_errors = false
                )
                if _timed_out(timeout, time_tmp[i], names[k], _abstols[k][j])
                    timed_out = true
                    break
                end
            end
            times[j, k] = timed_out ? NaN : mean(time_tmp) + x
            GC.gc()
        end
    end

    stats = nothing
    wps = [
        WorkPrecision(
            prob, _abstols[i], _reltols[i],
            _dicts_to_structarray(errors[i]),
            times[:, i], _dts[i], stats, names[i], error_estimate, N,
            _combined_tags(
                setups[i][:alg], _setup_tags(setups[i]),
                get(setups[i], :auto_tags, true)
            )
        )
            for i in 1:N
    ]
    return WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate,
        numruns_error,
        error_estimates === nothing ? Symbol[error_estimate] : collect(error_estimates)
    )
end

function WorkPrecisionSet(
        prob::AbstractEnsembleProblem, abstols, reltols, setups,
        test_dt = nothing;
        numruns = 5, trajectories = 1000,
        print_names = false, names = nothing, appxsol_setup = nothing,
        expected_value = nothing,
        error_estimate = :weak_final, error_estimates = nothing,
        ensemblealg = EnsembleThreads(), timeout = nothing,
        kwargs...
    )
    @assert names === nothing || length(setups) == length(names)

    weak_timeseries_errors = error_estimate ∈ WEAK_TIMESERIES_ERRORS
    weak_dense_errors = error_estimate ∈ WEAK_DENSE_ERRORS

    N = length(setups)
    M = length(abstols)
    times = Array{Float64}(undef, M, N)
    solutions = Array{Any}(undef, M, N)
    if names === nothing
        names = [_setup_name(setup) for setup in setups]
    end
    time_tmp = Vector{Float64}(undef, numruns)

    # First calculate all of the errors
    _abstols = [get(setups[k], :abstols, abstols) for k in 1:N]
    _reltols = [get(setups[k], :reltols, reltols) for k in 1:N]
    _dts = [get(setups[k], :dts, zeros(length(_abstols[k]))) for k in 1:N]
    for k in 1:N
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[k])

        for j in 1:M
            sol = solve(
                prob, setups[k][:alg], ensemblealg;
                filtered_setup...,
                abstol = _abstols[k][j],
                reltol = _reltols[k][j],
                dt = _dts[k][j],
                timeseries_errors = false,
                dense_errors = false,
                trajectories = Int(trajectories), kwargs...
            )
            solutions[j, k] = sol
        end
        @info "$(setups[k][:alg]) ($k/$N)"
    end

    if error_estimate ∈ WEAK_ERRORS
        if expected_value != nothing
            if error_estimate == :weak_final
                errors = [
                    [
                        LinearAlgebra.norm(
                            Statistics.mean(
                                solutions[i, j].u .-
                                    expected_value
                            )
                        )
                            for i in 1:M
                    ] for j in 1:N
                ]
            elseif error_estimate == :weak_l2
                errors = [
                    [
                        LinearAlgebra.norm(
                            Statistics.mean(
                                solutions[i, j] .-
                                    expected_value
                            )
                        )
                            for i in 1:M
                    ] for j in 1:N
                ]
            else
                error("Error estimate $error_estimate is not implemented yet.")
            end
        else
            sol = solve(
                prob, appxsol_setup[:alg], ensemblealg; kwargs..., appxsol_setup...,
                timeseries_errors = false, dense_errors = false,
                trajectories = Int(trajectories)
            )
            errors = [
                [
                    LinearAlgebra.norm(Statistics.mean(solutions[i, j].u .- sol.u))
                        for i in 1:M
                ] for j in 1:N
            ]
        end
    else
        error("use RODEProblem instead of EnsembleProblem for strong errors.")
    end

    local _sol

    # Now time it
    for k in 1:N
        # precompile
        GC.gc()
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[k])

        _sol = solve(
            prob, setups[k][:alg], ensemblealg;
            filtered_setup...,
            abstol = _abstols[k][1],
            reltol = _reltols[k][1],
            dt = _dts[k][1],
            timeseries_errors = false,
            dense_errors = false,
            trajectories = Int(trajectories), kwargs...
        )
        #x = isempty(_sol.t) ? 0 : round(Int,mean(_sol.t) - sum(_sol.t)/length(_sol.t))
        GC.gc()
        for j in 1:M
            timed_out = false
            for i in 1:numruns
                time_tmp[i] = @elapsed sol = solve(
                    prob, setups[k][:alg], ensemblealg;
                    filtered_setup...,
                    abstol = _abstols[k][j],
                    reltol = _reltols[k][j],
                    dt = _dts[k][j],
                    timeseries_errors = false,
                    dense_errors = false,
                    trajectories = Int(trajectories),
                    kwargs...
                )
                if _timed_out(timeout, time_tmp[i], names[k], _abstols[k][j])
                    timed_out = true
                    break
                end
            end
            times[j, k] = timed_out ? NaN : mean(time_tmp) #+ x
            GC.gc()
        end
    end
    stats = nothing
    wps = [
        WorkPrecision(
            prob, _abstols[i], _reltols[i],
            _dicts_to_structarray([Dict(error_estimate => err) for err in errors[i]]),
            times[:, i],
            _dts[i], stats, names[i], error_estimate, N,
            _combined_tags(
                setups[i][:alg], _setup_tags(setups[i]),
                get(setups[i], :auto_tags, true)
            )
        )
            for i in 1:N
    ]
    return WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate,
        Int(trajectories),
        error_estimates === nothing ? Symbol[error_estimate] : collect(error_estimates)
    )
end

function WorkPrecisionSet(
        prob::AbstractBVProblem,
        abstols, reltols, setups;
        print_names = false, names = nothing, appxsol = nothing,
        error_estimate = :final, error_estimates = nothing,
        test_dt = nothing, timeout = nothing, kwargs...
    )
    N = length(setups)
    @assert names === nothing || length(setups) == length(names)
    wps = Vector{WorkPrecision}(undef, N)
    if names === nothing
        names = [_setup_name(setup) for setup in setups]
    end

    active = error_estimates === nothing ? Symbol[error_estimate] : collect(error_estimates)
    timeseries_errors = error_estimates === nothing ? nothing : _needs_timeseries(active)
    dense_errors = error_estimates === nothing ? nothing : _needs_dense(active)

    for i in 1:N
        print_names && println(names[i])
        _abstols = get(setups[i], :abstols, abstols)
        _reltols = get(setups[i], :reltols, reltols)
        _dts = get(setups[i], :dts, nothing)
        filtered_setup = filter(p -> p.first in SciMLBase.allowedkeywords, setups[i])

        wps[i] = WorkPrecision(
            prob, setups[i][:alg], _abstols, _reltols, _dts;
            appxsol,
            error_estimate,
            timeout,
            timeseries_errors,
            dense_errors,
            tags = _setup_tags(setups[i]),
            auto_tags = get(setups[i], :auto_tags, true),
            name = names[i], kwargs..., filtered_setup...
        )
    end
    return WorkPrecisionSet(
        wps, N, abstols, reltols, prob, setups, names, error_estimate,
        nothing, active
    )
end

"""
    get_sample_errors(
        prob::AbstractRODEProblem, setup, test_dt = nothing;
        numruns, solution_runs, appxsol_setup = nothing,
        sample_error_runs = 10^7, parallel_type = :none, kwargs...
    )

Estimate an approximate 95% confidence half-width for the sampling error of an RODE
solver setup. `numruns` may be one sample count or a collection of counts; the return
value is respectively a scalar or a collection of sampling-error estimates.
`solution_runs` controls the number of independent estimates used to measure their
variation.

When an analytic solution is available, `sample_error_runs` controls the Monte Carlo
estimate of its expected endpoint. Otherwise, repeated numerical solution means provide
the endpoint reference. Set `parallel_type = :threads` to parallelize the solution
samples.
"""
function get_sample_errors(
        prob::AbstractRODEProblem, setup, test_dt = nothing;
        appxsol_setup = nothing,
        numruns, error_estimate = :final,
        sample_error_runs = Int(1.0e7),
        solution_runs,
        parallel_type = :none, kwargs...
    )
    maxnumruns = findmax(numruns)[1]

    tmp_solutions_full = map(1:solution_runs) do i
        @info "Solution Run: $i"
        # Use the WorkPrecision stuff to calculate the errors
        tmp_solutions = Array{Any}(undef, maxnumruns, 1, 1)
        setups = [setup]
        abstols = [1.0e-2] # Standard default
        reltols = [1.0e-2] # Standard default
        M = 1
        N = 1
        timeseries_errors = false
        dense_errors = false
        if parallel_type == :threads
            Threads.@threads for i in 1:maxnumruns
                @error_calculation
            end
        elseif parallel_type == :none
            for i in 1:maxnumruns
                @error_calculation
            end
        end
        tmp_solutions = vec(tmp_solutions)
    end

    if SciMLBase.has_analytic(prob.f)
        analytical_mean_end = mean(1:sample_error_runs) do i
            _dt = prob.tspan[2] - prob.tspan[1]
            if prob.u0 isa Number
                W = sqrt(_dt) * randn()
            else
                W = sqrt(_dt) * randn(size(prob.u0))
            end
            prob.f.analytic(prob.u0, prob.p, prob.tspan[2], W)
        end
    else
        # Use the mean of the means as the analytical mean
        analytical_mean_end = mean(
            mean(
                tmp_solutions[i].u[end]
                    for i in 1:length(tmp_solutions)
            )
                for tmp_solutions in tmp_solutions_full
        )
    end

    if numruns isa Number
        mean_solution_ends = [
            mean([tmp_solutions[i].u[end] for i in 1:maxnumruns])
                for tmp_solutions in tmp_solutions_full
        ]
        return sample_error = 1.96std(
            norm(mean_sol_end - analytical_mean_end)
                for mean_sol_end in mean_solution_ends
        ) /
            sqrt(numruns)
    else
        map(1:length(numruns)) do i
            mean_solution_ends = [
                mean([tmp_solutions[i].u[end] for i in 1:numruns[i]])
                    for tmp_solutions in tmp_solutions_full
            ]
            sample_error = 1.96std(
                norm(mean_sol_end - analytical_mean_end)
                    for mean_sol_end in mean_solution_ends
            ) /
                sqrt(numruns[i])
        end
    end
end

## Tagging and filtering

_as_tags(tags) = tags isa Symbol ? (tags,) : Tuple(tags)

"""
    get_tags(wp_set::WorkPrecisionSet) -> Vector{Vector{Symbol}}

Return the tags of each entry of `wp_set`, in order. Tags come from the `:tags` entry
of the corresponding setup dictionary.
"""
get_tags(wp_set::WorkPrecisionSet) = [wp.tags for wp in wp_set.wps]

"""
    unique_tags(wp_set::WorkPrecisionSet) -> Vector{Symbol}

Return every tag used by any entry of `wp_set`, sorted and deduplicated.
"""
function unique_tags(wp_set::WorkPrecisionSet)
    alltags = Symbol[]
    for wp in wp_set.wps
        append!(alltags, wp.tags)
    end
    return unique!(sort!(alltags))
end

_has_all_tags(wp, tags) = all(t -> t in wp.tags, tags)
_has_any_tag(wp, tags) = any(t -> t in wp.tags, tags)

"""
    filter_by_tags(wp_set::WorkPrecisionSet, tags::Symbol...) -> WorkPrecisionSet

Return the entries of `wp_set` tagged with every one of `tags` (AND logic), as a new
`WorkPrecisionSet`. Passing no tags returns `wp_set` unchanged.
"""
function filter_by_tags(wp_set::WorkPrecisionSet, tags::Symbol...)
    isempty(tags) && return wp_set
    return _subset_wps(wp_set, findall(wp -> _has_all_tags(wp, tags), wp_set.wps))
end

"""
    exclude_by_tags(wp_set::WorkPrecisionSet, tags::Symbol...) -> WorkPrecisionSet

Return the entries of `wp_set` carrying none of `tags` (OR logic on the exclusion), as
a new `WorkPrecisionSet`. Passing no tags returns `wp_set` unchanged.
"""
function exclude_by_tags(wp_set::WorkPrecisionSet, tags::Symbol...)
    isempty(tags) && return wp_set
    return _subset_wps(wp_set, findall(wp -> !_has_any_tag(wp, tags), wp_set.wps))
end

"""
    merge_wp_sets(sets::WorkPrecisionSet...) -> WorkPrecisionSet

Concatenate the entries of several `WorkPrecisionSet`s into one. The tolerances,
problem and error estimates are taken from the first set, so merging results computed
on different problems or tolerance grids gives a set whose metadata describes only the
first of them.
"""
function merge_wp_sets(sets::WorkPrecisionSet...)
    isempty(sets) && throw(ArgumentError("At least one WorkPrecisionSet is required"))
    wps = reduce(vcat, [s.wps for s in sets])
    setups = reduce(vcat, [s.setups for s in sets])
    names = reduce(vcat, [s.names for s in sets])
    first_set = first(sets)
    return WorkPrecisionSet(
        wps, length(wps), first_set.abstols, first_set.reltols, first_set.prob,
        setups, names, first_set.error_estimate, first_set.numruns,
        first_set.active_error_estimates
    )
end

function _subset_wps(wp_set::WorkPrecisionSet, indices)
    isempty(indices) &&
        @warn "No entries match the requested tags. Returning an empty WorkPrecisionSet."
    setups = wp_set.setups isa AbstractVector ? wp_set.setups[indices] : wp_set.setups
    names = wp_set.names isa AbstractVector ? wp_set.names[indices] : wp_set.names
    return WorkPrecisionSet(
        wp_set.wps[indices], length(indices), wp_set.abstols, wp_set.reltols,
        wp_set.prob, setups, names, wp_set.error_estimate, wp_set.numruns,
        wp_set.active_error_estimates
    )
end

# Indices selected by the plot recipe's tag keywords: `tags` narrows (AND), then
# `include_tags` adds back entries matching all of those tags, then `exclude_tags`
# drops any entry carrying one of them.
function _selected_indices(
        wp_set::WorkPrecisionSet; tags = nothing, include_tags = nothing,
        exclude_tags = nothing
    )
    indices = collect(eachindex(wp_set.wps))
    if tags !== nothing
        ts = _as_tags(tags)
        indices = filter(i -> _has_all_tags(wp_set.wps[i], ts), indices)
    end
    if include_tags !== nothing
        ts = _as_tags(include_tags)
        extra = findall(wp -> _has_all_tags(wp, ts), wp_set.wps)
        indices = sort!(union(indices, extra))
    end
    if exclude_tags !== nothing
        ts = _as_tags(exclude_tags)
        indices = filter(i -> !_has_any_tag(wp_set.wps[i], ts), indices)
    end
    return indices
end

## Best-of-family selection

# (log10(error), log10(time)) for the tolerances that produced a usable measurement,
# sorted by increasing error.
function _wp_points(wp::WorkPrecision)
    errs = getproperty(wp.errors, wp.error_estimate)
    pts = [
        (log10(e), log10(t)) for (e, t) in zip(errs, wp.times)
            if !isnan(e) && !isnan(t) && e > 0 && t > 0
    ]
    return sort!(pts; by = first)
end

function _trapezoid(pts)
    return sum(
        0.5 * (pts[i][2] + pts[i - 1][2]) * (pts[i][1] - pts[i - 1][1])
            for i in 2:length(pts)
    )
end

"""
    wp_area(wp::WorkPrecision) -> Float64

Trapezoidal area under the log₁₀(time)-vs-log₁₀(error) curve of `wp`, using the
tolerances that produced a usable (positive, non-`NaN`) error and time. Lower is
better. Returns `Inf` when fewer than two such points exist.

The area grows with the width of the error range covered, so it only compares methods
that span a similar range; [`best_by_tag`](@ref) normalizes by that width instead.
"""
function wp_area(wp::WorkPrecision)
    pts = _wp_points(wp)
    length(pts) < 2 && return Inf
    return _trapezoid(pts)
end

# Rank key: most usable tolerance points first, then lowest mean log10 time across the
# error range covered. Ranking on the mean rather than the raw area keeps methods that
# reach a wider range of accuracies from being penalized for it.
function _wp_rank(wp::WorkPrecision)
    pts = _wp_points(wp)
    length(pts) < 2 && return (-length(pts), Inf)
    span = pts[end][1] - pts[1][1]
    span == 0 && return (-length(pts), mean(p[2] for p in pts))
    return (-length(pts), _trapezoid(pts) / span)
end

function _best_indices(wp_set::WorkPrecisionSet, tag::Symbol, n::Int, metric::Symbol)
    metric === :area ||
        throw(ArgumentError("Unknown metric: $metric. Supported metrics: :area"))
    candidates = findall(wp -> tag in wp.tags, wp_set.wps)
    isempty(candidates) && return candidates
    perm = sortperm([_wp_rank(wp_set.wps[i]) for i in candidates])
    return candidates[perm[1:min(n, length(perm))]]
end

"""
    best_by_tag(wp_set::WorkPrecisionSet, tag::Symbol; n = 1, metric = :area)
        -> WorkPrecisionSet

Return the `n` best-performing entries of `wp_set` tagged `tag`. With `metric = :area`
(the only metric currently supported) entries are ranked first by how many tolerances
produced a usable measurement — so a method that fails at tight tolerances does not win
on the few points it survived — and then by their mean log₁₀ solve time over the
log₁₀ error range they cover, i.e. [`wp_area`](@ref) normalized by that range.
"""
function best_by_tag(
        wp_set::WorkPrecisionSet, tag::Symbol; n::Int = 1, metric::Symbol = :area
    )
    return _subset_wps(wp_set, _best_indices(wp_set, tag, n, metric))
end

"""
    best_of_families(wp_set::WorkPrecisionSet, family_tags; n = 1, metric = :area)
        -> WorkPrecisionSet

Combine the `n` best entries of each family in `family_tags` (see
[`best_by_tag`](@ref)) into one `WorkPrecisionSet` for a cross-family comparison. An
entry belonging to several families is included once. Throws an `ArgumentError` when no
entry carries any of `family_tags`.
"""
function best_of_families(
        wp_set::WorkPrecisionSet, family_tags; n::Int = 1, metric::Symbol = :area
    )
    indices = _best_family_indices(wp_set, family_tags, n, metric)
    isempty(indices) &&
        throw(ArgumentError("No entries found for any of the family tags $family_tags"))
    return _subset_wps(wp_set, indices)
end

function _best_family_indices(wp_set::WorkPrecisionSet, family_tags, n, metric)
    indices = Int[]
    for tag in family_tags
        append!(indices, _best_indices(wp_set, tag, n, metric))
    end
    return sort!(unique!(indices))
end

## AutoDiff comparison helpers

# Short name of an AD backend, e.g. AutoForwardDiff() -> "ForwardDiff".
function ad_backend_name(backend)
    name = string(nameof(typeof(backend)))
    return startswith(name, "Auto") ? name[5:end] : name
end

"""
    with_autodiff_variants(setups; ad_backends, tag_prefix = :autodiff)
        -> Vector{Dict{Symbol, Any}}

Expand `setups` with one copy per entry of `ad_backends`, each copy replacing `:alg`
with the same algorithm rebuilt for that backend. Setups whose algorithm takes no
`autodiff` argument (explicit Runge-Kutta methods, say) are passed through unexpanded.

The originals are kept and tagged `\$(tag_prefix)_default`; each variant is tagged with
the lowercased backend name, e.g. `:autodiff_forwarddiff` for `AutoForwardDiff()`, and
named after its backend so the legends stay readable. Existing tags are preserved and
the input setups are not mutated, so the result plots as an AD comparison via
`plot(wp_set, tags = [:autodiff_forwarddiff])`.
"""
function with_autodiff_variants(setups; ad_backends, tag_prefix::Symbol = :autodiff)
    result = Vector{Dict{Symbol, Any}}()
    for setup in setups
        original = Dict{Symbol, Any}(setup)
        original[:tags] = [_setup_tags(setup); Symbol(tag_prefix, :_default)]
        push!(result, original)

        hasfield(typeof(setup[:alg]), :autodiff) || continue
        for backend in ad_backends
            variant = Dict{Symbol, Any}(setup)
            variant[:alg] = remake(setup[:alg]; autodiff = backend)
            variant[:name] = "$(_setup_name(setup)) ($(ad_backend_name(backend)))"
            variant[:tags] = [
                _setup_tags(setup);
                Symbol(tag_prefix, :_, lowercase(ad_backend_name(backend)))
            ]
            push!(result, variant)
        end
    end
    return result
end

Base.length(wp::WorkPrecision) = wp.N
Base.size(wp::WorkPrecision) = length(wp)
Base.getindex(wp::WorkPrecision, i::Int) = wp.times[i]
Base.getindex(wp::WorkPrecision, ::Colon) = wp.times
Base.firstindex(wp::WorkPrecision) = 1
Base.lastindex(wp::WorkPrecision) = lastindex(wp.times)

function Base.show(io::IO, wp::WorkPrecision)
    println(io, "Name: $(wp.name)")
    println(io, "Times: $(wp.times)")
    return println(io, "Errors: $(wp.errors)")
end

Base.length(wp_set::WorkPrecisionSet) = wp_set.N
Base.size(wp_set::WorkPrecisionSet) = length(wp_set)
Base.getindex(wp_set::WorkPrecisionSet, i::Int) = wp_set.wps[i]
Base.getindex(wp_set::WorkPrecisionSet, ::Colon) = wp_set.wps
Base.firstindex(wp_set::WorkPrecisionSet) = 1
Base.lastindex(wp_set::WorkPrecisionSet) = lastindex(wp_set.wps)

function Base.show(io::IO, wp_set::WorkPrecisionSet)
    return println(io, "WorkPrecisionSet of $(wp_set.N) wps")
end

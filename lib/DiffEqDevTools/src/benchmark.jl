using Statistics

# Default names of algorithms:
# Workaround for `MethodOfSteps` algorithms, otherwise they are all called "MethodOfSteps"
# Ideally this would be a trait (in SciMLBase?), so packages could implement it
_default_name(alg) = __default_name(alg)
function _default_name(alg::AbstractDDEAlgorithm)
    isdefined(alg, :alg) ? _default_name(alg.alg) : __default_name(alg)
end
__default_name(alg) = string(nameof(typeof(alg)))

## Shootouts

mutable struct Shootout
    setups::Vector{Dict{Symbol, Any}}
    times::Any#::Vector{Float64}
    errors::Any#::Vector{uType}
    effs::Any#::Vector{Float64} # Efficiencies
    effratios::Any#::Matrix{uEltype}
    solutions::Any
    names::Vector{String}
    N::Int
    bestidx::Int
    winner::String
end

mutable struct ShootoutSet
    shootouts::Vector{Shootout}
    probs::Any#::Vector{DEProblem}
    probaux::Any#::Vector{Dict{Symbol,Any}}
    N::Int
    winners::Vector{String}
end

function ode_shootout(args...; kwargs...)
    @warn("ode_shootout is deprecated. Use ShootOut instead")
    ShootOut(args...; kwargs...)
end

function Shootout(
        prob, setups; appxsol = nothing, names = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, kwargs...)
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
        sol = solve(prob, setups[i][:alg]; timeseries_errors = timeseries_errors,
            dense_errors = dense_errors, kwargs..., setups[i]...)

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
            () -> @elapsed solve(_prob, alg, sol.u, sol.t, sol.k;
                timeseries_errors = false,
                dense_errors = false, kwargs...)
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
    return Shootout(setups, times, errors, effs, effratios, solutions, names, N, bestidx,
        winner)
end

function ode_shootoutset(args...; kwargs...)
    @warn("ode_shootoutset is deprecated. Use ShootoutSet instead")
    ShootoutSet(args...; kwargs...)
end

function ShootoutSet(probs, setups; probaux = nothing,
        names = nothing, print_names = false, kwargs...)
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
        shootouts[i] = Shootout(probs[i], setups; names = names, kwargs..., probaux[i]...)
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
    println(io, "EffRatios: $(shoot.effratios[shoot.bestidx,:])")
end

Base.length(set::ShootoutSet) = set.N
Base.size(set::ShootoutSet) = length(set)
Base.getindex(set::ShootoutSet, i::Int) = set.shootouts[i]
Base.getindex(set::ShootoutSet, ::Colon) = set.shootouts
Base.show(io::IO, set::ShootoutSet) = print(io, "ShootoutSet of $(set.N) shootouts ")
Base.firstindex(shoot::ShootoutSet) = 1
Base.lastindex(shoot::ShootoutSet) = lastindex(shoot.shootouts)

## WorkPrecisions

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
end

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
end

function WorkPrecision(prob, alg, abstols, reltols, dts = nothing;
        name = nothing, appxsol = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, kwargs...)
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
        timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
        dense_errors = error_estimate ∈ DENSE_ERRORS
        for i in 1:N
            if dts === nothing
                sol = solve(_prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
            else
                sol = solve(_prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], dt = dts[i],
                    timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
            end

            stats[i] = sol.stats

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
                            () -> @elapsed solve(_prob, alg, sol.u, sol.t;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        else
                            () -> @elapsed solve(_prob, alg, sol.u, sol.t, sol.k;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        end
                    else
                        if _prob isa DAEProblem
                            () -> @elapsed solve(_prob, alg, sol.u, sol.t;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        else
                            () -> @elapsed solve(_prob, alg, sol.u, sol.t, sol.k;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
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
                errors[i] = Dict(
                    :l∞ => NaN, :L2 => NaN, :final => NaN, :l2 => NaN, :L∞ => NaN)
                times[i] = NaN
            end
        end
    end
    return WorkPrecision(prob, abstols, reltols, StructArray(NamedTuple.(errors)),
        times, dts, stats, name, error_estimate, N)
end

# Work precision information for a BVP
function WorkPrecision(prob::AbstractBVProblem, alg, abstols, reltols, dts = nothing;
        name = nothing, appxsol = nothing, error_estimate = :final,
        numruns = 20, seconds = 2, kwargs...)
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
        timeseries_errors = error_estimate ∈ TIMESERIES_ERRORS
        dense_errors = error_estimate ∈ DENSE_ERRORS
        for i in 1:N
            if dts === nothing
                sol = solve(_prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
            else
                sol = solve(_prob, alg; kwargs..., abstol = abstols[i],
                    reltol = reltols[i], dt = dts[i],
                    timeseries_errors = timeseries_errors,
                    dense_errors = dense_errors)
            end

            stats[i] = sol.stats
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
                            () -> @elapsed solve(_prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        else
                            () -> @elapsed solve(_prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        end
                    else
                        if _prob isa DAEProblem
                            () -> @elapsed solve(_prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
                        else
                            () -> @elapsed solve(_prob, alg;
                                abstol = abstols[i],
                                reltol = reltols[i],
                                dt = dts[i],
                                timeseries_errors = false,
                                dense_errors = false, kwargs...)
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
                errors[i] = Dict(
                    :l∞ => NaN, :L2 => NaN, :final => NaN, :l2 => NaN, :L∞ => NaN)
                times[i] = NaN
            end
        end
    end
    return WorkPrecision(prob, abstols, reltols, StructArray(NamedTuple.(errors)),
        times, dts, stats, name, error_estimate, N)
end

# Work precision information for a nonlinear problem.
function WorkPrecision(
        prob::NonlinearProblem, alg, abstols, reltols, dts = nothing; name = nothing,
        appxsol = nothing, error_estimate = :l2, numruns = 20, seconds = 2, kwargs...)
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
            sol = solve(_prob, alg; kwargs..., abstol = abstols[i], reltol = reltols[i])

            stats[i] = sol.stats

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

                () -> @elapsed solve(_prob, alg;
                    abstol = abstols[i],
                    reltol = reltols[i],
                    kwargs...)
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

    return WorkPrecision(prob, abstols, reltols, StructArray(NamedTuple.(errors)),
        times, dts, stats, name, error_estimate, N)
end

function WorkPrecisionSet(prob,
        abstols, reltols, setups;
        print_names = false, names = nothing, appxsol = nothing,
        error_estimate = :final,
        test_dt = nothing, kwargs...)
    N = length(setups)
    @assert names === nothing || length(setups) == length(names)
    wps = Vector{WorkPrecision}(undef, N)
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
    end
    for i in 1:N
        print_names && println(names[i])
        _abstols = get(setups[i], :abstols, abstols)
        _reltols = get(setups[i], :reltols, reltols)
        _dts = get(setups[i], :dts, nothing)
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[i])

        wps[i] = WorkPrecision(prob, setups[i][:alg], _abstols, _reltols, _dts;
            appxsol = appxsol,
            error_estimate = error_estimate,
            name = names[i], kwargs..., filtered_setup...)
    end
    return WorkPrecisionSet(wps, N, abstols, reltols, prob, setups, names, error_estimate,
        nothing)
end

@def error_calculation begin
    if !DiffEqBase.has_analytic(prob.f)
        t = prob.tspan[1]:test_dt:prob.tspan[2]
        brownian_values = cumsum([[zeros(size(prob.u0))];
                                  [sqrt(test_dt) * randn(size(prob.u0))
                                   for i in 1:(length(t) - 1)]])
        brownian_values2 = cumsum([[zeros(size(prob.u0))];
                                   [sqrt(test_dt) * randn(size(prob.u0))
                                    for i in 1:(length(t) - 1)]])
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
    filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[1])

    sol = solve(_prob, setups[1][:alg];
        kwargs..., filtered_setup..., abstol = _abstols[1],
        reltol = _reltols[1], dt = _dts[1],
        timeseries_errors = false,
        dense_errors = false)

    for j in 1:M, k in 1:N
        _abstols = get(setups[k], :abstols, abstols)
        _reltols = get(setups[k], :reltols, reltols)
        _dts = get(setups[k], :dts, zeros(length(_abstols)))
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[k])

        sol = solve(_prob, setups[k][:alg];
            kwargs..., filtered_setup..., abstol = _abstols[j],
            reltol = _reltols[j], dt = _dts[j],
            timeseries_errors = timeseries_errors,
            dense_errors = dense_errors)
        DiffEqBase.has_analytic(prob.f) ? err_sol = sol : err_sol = appxtrue(sol, true_sol)
        tmp_solutions[i, j, k] = err_sol
    end
end

function WorkPrecisionSet(prob::AbstractRODEProblem, abstols, reltols, setups,
        test_dt = nothing;
        numruns = 20, numruns_error = 20,
        print_names = false, names = nothing, appxsol_setup = nothing,
        error_estimate = :final, parallel_type = :none,
        kwargs...)
    @assert names === nothing || length(setups) == length(names)
    timeseries_errors = DiffEqBase.has_analytic(prob.f) &&
                        error_estimate ∈ TIMESERIES_ERRORS
    weak_timeseries_errors = error_estimate ∈ WEAK_TIMESERIES_ERRORS
    weak_dense_errors = error_estimate ∈ WEAK_DENSE_ERRORS
    dense_errors = DiffEqBase.has_analytic(prob.f) && error_estimate ∈ DENSE_ERRORS
    N = length(setups)
    M = length(abstols)
    times = Array{Float64}(undef, M, N)
    tmp_solutions = Array{Any}(undef, numruns_error, M, N)
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
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

    _solutions_k = [[EnsembleSolution(tmp_solutions[:, j, k], 0.0, true) for j in 1:M]
                    for k in 1:N]
    solutions = [[DiffEqBase.calculate_ensemble_errors(sim;
                      weak_timeseries_errors = weak_timeseries_errors,
                      weak_dense_errors = weak_dense_errors)
                  for sim in sol_k] for sol_k in _solutions_k]
    if error_estimate ∈ WEAK_ERRORS
        errors = [[solutions[j][i].weak_errors for i in 1:M] for j in 1:N]
    else
        errors = [[solutions[j][i].error_means for i in 1:M] for j in 1:N]
    end

    local _sol

    # Now time it
    _abstols = [get(setups[k], :abstols, abstols) for k in 1:N]
    _reltols = [get(setups[k], :reltols, reltols) for k in 1:N]
    _dts = [get(setups[k], :dts, zeros(length(_abstols))) for k in 1:N]
    for k in 1:N
        # precompile
        GC.gc()
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[k])

        _sol = solve(prob, setups[k][:alg];
            kwargs..., filtered_setup..., abstol = _abstols[k][1],
            reltol = _reltols[k][1], dt = _dts[k][1],
            timeseries_errors = false,
            dense_errors = false)
        x = isempty(_sol.t) ? 0 : round(Int, mean(_sol.t) - sum(_sol.t) / length(_sol.t))
        GC.gc()
        for j in 1:M
            for i in 1:numruns
                time_tmp[i] = @elapsed sol = solve(prob, setups[k][:alg];
                    kwargs..., filtered_setup...,
                    abstol = _abstols[k][j],
                    reltol = _reltols[k][j], dt = _dts[k][j],
                    timeseries_errors = false,
                    dense_errors = false)
            end
            times[j, k] = mean(time_tmp) + x
            GC.gc()
        end
    end

    stats = nothing
    wps = [WorkPrecision(prob, _abstols[i], _reltols[i],
               StructArray(NamedTuple.(errors[i])),
               times[:, i], _dts[i], stats, names[i], error_estimate, N)
           for i in 1:N]
    WorkPrecisionSet(wps, N, abstols, reltols, prob, setups, names, error_estimate,
        numruns_error)
end

function WorkPrecisionSet(prob::AbstractEnsembleProblem, abstols, reltols, setups,
        test_dt = nothing;
        numruns = 5, trajectories = 1000,
        print_names = false, names = nothing, appxsol_setup = nothing,
        expected_value = nothing,
        error_estimate = :weak_final, ensemblealg = EnsembleThreads(),
        kwargs...)
    @assert names === nothing || length(setups) == length(names)

    weak_timeseries_errors = error_estimate ∈ WEAK_TIMESERIES_ERRORS
    weak_dense_errors = error_estimate ∈ WEAK_DENSE_ERRORS

    N = length(setups)
    M = length(abstols)
    times = Array{Float64}(undef, M, N)
    solutions = Array{Any}(undef, M, N)
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
    end
    time_tmp = Vector{Float64}(undef, numruns)

    # First calculate all of the errors
    _abstols = [get(setups[k], :abstols, abstols) for k in 1:N]
    _reltols = [get(setups[k], :reltols, reltols) for k in 1:N]
    _dts = [get(setups[k], :dts, zeros(length(_abstols))) for k in 1:N]
    for k in 1:N
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[k])

        for j in 1:M
            sol = solve(prob, setups[k][:alg], ensemblealg;
                filtered_setup...,
                abstol = _abstols[k][j],
                reltol = _reltols[k][j],
                dt = _dts[k][j],
                timeseries_errors = false,
                dense_errors = false,
                trajectories = Int(trajectories), kwargs...)
            solutions[j, k] = sol
        end
        @info "$(setups[k][:alg]) ($k/$N)"
    end

    if error_estimate ∈ WEAK_ERRORS
        if expected_value != nothing
            if error_estimate == :weak_final
                errors = [[LinearAlgebra.norm(Statistics.mean(solutions[i, j].u .-
                                                              expected_value))
                           for i in 1:M] for j in 1:N]
            elseif error_estimate == :weak_l2
                errors = [[LinearAlgebra.norm(Statistics.mean(solutions[i, j] .-
                                                              expected_value))
                           for i in 1:M] for j in 1:N]
            else
                error("Error estimate $error_estimate is not implemented yet.")
            end
        else
            sol = solve(
                prob, appxsol_setup[:alg], ensemblealg; kwargs..., appxsol_setup...,
                timeseries_errors = false, dense_errors = false,
                trajectories = Int(trajectories))
            errors = [[LinearAlgebra.norm(Statistics.mean(solutions[i, j].u .- sol.u))
                       for i in 1:M] for j in 1:N]
        end
    else
        error("use RODEProblem instead of EnsembleProblem for strong errors.")
    end

    local _sol

    # Now time it
    for k in 1:N
        # precompile
        GC.gc()
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[k])

        _sol = solve(prob, setups[k][:alg], ensemblealg;
            filtered_setup...,
            abstol = _abstols[k][1],
            reltol = _reltols[k][1],
            dt = _dts[k][1],
            timeseries_errors = false,
            dense_errors = false,
            trajectories = Int(trajectories), kwargs...)
        #x = isempty(_sol.t) ? 0 : round(Int,mean(_sol.t) - sum(_sol.t)/length(_sol.t))
        GC.gc()
        for j in 1:M
            for i in 1:numruns
                time_tmp[i] = @elapsed sol = solve(prob, setups[k][:alg], ensemblealg;
                    filtered_setup...,
                    abstol = _abstols[k][j],
                    reltol = _reltols[k][j],
                    dt = _dts[k][j],
                    timeseries_errors = false,
                    dense_errors = false,
                    trajectories = Int(trajectories),
                    kwargs...)
            end
            times[j, k] = mean(time_tmp) #+ x
            GC.gc()
        end
    end
    stats = nothing
    wps = [WorkPrecision(prob, _abstols[i], _reltols[i], errors[i], times[:, i],
               _dts[i], stats, names[i], error_estimate, N)
           for i in 1:N]
    WorkPrecisionSet(wps, N, abstols, reltols, prob, setups, names, error_estimate,
        Int(trajectories))
end

function WorkPrecisionSet(prob::AbstractBVProblem,
        abstols, reltols, setups;
        print_names = false, names = nothing, appxsol = nothing,
        error_estimate = :final,
        test_dt = nothing, kwargs...)
    N = length(setups)
    @assert names === nothing || length(setups) == length(names)
    wps = Vector{WorkPrecision}(undef, N)
    if names === nothing
        names = [_default_name(setup[:alg]) for setup in setups]
    end
    for i in 1:N
        print_names && println(names[i])
        _abstols = get(setups[i], :abstols, abstols)
        _reltols = get(setups[i], :reltols, reltols)
        _dts = get(setups[i], :dts, nothing)
        filtered_setup = filter(p -> p.first in DiffEqBase.allowedkeywords, setups[i])

        wps[i] = WorkPrecision(prob, setups[i][:alg], _abstols, _reltols, _dts;
            appxsol = appxsol,
            error_estimate = error_estimate,
            name = names[i], kwargs..., filtered_setup...)
    end
    return WorkPrecisionSet(wps, N, abstols, reltols, prob, setups, names, error_estimate,
        nothing)
end

function get_sample_errors(prob::AbstractRODEProblem, setup, test_dt = nothing;
        appxsol_setup = nothing,
        numruns, error_estimate = :final,
        sample_error_runs = Int(1e7),
        solution_runs,
        parallel_type = :none, kwargs...)
    maxnumruns = findmax(numruns)[1]

    tmp_solutions_full = map(1:solution_runs) do i
        @info "Solution Run: $i"
        # Use the WorkPrecision stuff to calculate the errors
        tmp_solutions = Array{Any}(undef, maxnumruns, 1, 1)
        setups = [setup]
        abstols = [1e-2] # Standard default
        reltols = [1e-2] # Standard default
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

    if DiffEqBase.has_analytic(prob.f)
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
        analytical_mean_end = mean(mean(tmp_solutions[i].u[end]
                                   for i in 1:length(tmp_solutions))
        for tmp_solutions in tmp_solutions_full)
    end

    if numruns isa Number
        mean_solution_ends = [mean([tmp_solutions[i].u[end] for i in 1:maxnumruns])
                              for tmp_solutions in tmp_solutions_full]
        return sample_error = 1.96std(norm(mean_sol_end - analytical_mean_end)
        for mean_sol_end in mean_solution_ends) /
                              sqrt(numruns)
    else
        map(1:length(numruns)) do i
            mean_solution_ends = [mean([tmp_solutions[i].u[end] for i in 1:numruns[i]])
                                  for tmp_solutions in tmp_solutions_full]
            sample_error = 1.96std(norm(mean_sol_end - analytical_mean_end)
            for mean_sol_end in mean_solution_ends) /
                           sqrt(numruns[i])
        end
    end
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
    println(io, "Errors: $(wp.errors)")
end

Base.length(wp_set::WorkPrecisionSet) = wp_set.N
Base.size(wp_set::WorkPrecisionSet) = length(wp_set)
Base.getindex(wp_set::WorkPrecisionSet, i::Int) = wp_set.wps[i]
Base.getindex(wp_set::WorkPrecisionSet, ::Colon) = wp_set.wps
Base.firstindex(wp_set::WorkPrecisionSet) = 1
Base.lastindex(wp_set::WorkPrecisionSet) = lastindex(wp_set.wps)

function Base.show(io::IO, wp_set::WorkPrecisionSet)
    println(io, "WorkPrecisionSet of $(wp_set.N) wps")
end

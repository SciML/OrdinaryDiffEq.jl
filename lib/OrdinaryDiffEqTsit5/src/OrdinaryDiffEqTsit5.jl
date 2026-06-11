module OrdinaryDiffEqTsit5

import OrdinaryDiffEqCore: alg_stability_size, explicit_rk_docstring,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
    alg_cache,
    OrdinaryDiffEqConstantCache, @fold, trivial_limiter!,
    constvalue, perform_step!, @cache,
    _ode_interpolant, _ode_interpolant!,
    CompiledFloats, @OnDemandTableauExtract,
    CompositeAlgorithm, _ode_addsteps!,
    AutoAlgSwitch, get_fsalfirstlast,
    full_cache, DerivativeOrderNotPossibleError,
    TmpCache, build_tmp_cache
using FastBroadcast: Serial
import MuladdMacro: @muladd
import FastBroadcast: @..
import RecursiveArrayTools: recursivefill!, recursive_unitless_bottom_eltype,
    copyat_or_push!
import LinearAlgebra: norm
using TruncatedStacktraces: @truncate_stacktrace
import SciMLBase: alg_order, @def
using DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: initialize!
import OrdinaryDiffEqCore
using CommonSolve: solve

using Reexport: Reexport, @reexport
@reexport using SciMLBase
using SciMLBase: SciMLBase, ODEProblem

include("algorithms.jl")
include("alg_utils.jl")
include("tsit_caches.jl")
include("tsit_tableaus.jl")
include("interp_func.jl")
include("interpolants.jl")
include("tsit_perform_step.jl")

import PrecompileTools
import Preferences

PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [Tsit5()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileAutoDePSpecialize", true)
        push!(
            prob_list,
            ODEProblem{true, OrdinaryDiffEqCore.AutoDePSpecialize}(
                OrdinaryDiffEqCore.lorenz_p, [1.0; 0.0; 0.0],
                (0.0, 1.0), OrdinaryDiffEqCore.lorenz_p_params
            )
        )
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(
                lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]
            )
        )
    end

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export Tsit5, AutoTsit5

# Cross-sublibrary cache types that other OrdinaryDiffEq solver sublibraries
# (e.g. OrdinaryDiffEqNordsieck) reference to reuse the Tsit5 step. Marked
# public so those references are recognized as a supported extension API rather
# than internal access.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Tsit5Cache, :Tsit5ConstantCache))
end

end

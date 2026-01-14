module OrdinaryDiffEqTsit5

import OrdinaryDiffEqCore: ODEVerbosity, alg_order, alg_stability_size, explicit_rk_docstring,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
    alg_cache,
    OrdinaryDiffEqConstantCache, @fold, trivial_limiter!,
    constvalue, perform_step!, calculate_residuals, @cache,
    calculate_residuals!, _ode_interpolant, _ode_interpolant!,
    CompiledFloats, @OnDemandTableauExtract, initialize!,
    perform_step!,
    CompositeAlgorithm, _ode_addsteps!, copyat_or_push!,
    AutoAlgSwitch, get_fsalfirstlast,
    full_cache, DerivativeOrderNotPossibleError
import Static: False
import MuladdMacro: @muladd
import FastBroadcast: @..
import RecursiveArrayTools: recursivefill!, recursive_unitless_bottom_eltype
import LinearAlgebra: norm
using TruncatedStacktraces: @truncate_stacktrace
import SciMLBase: @def
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

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

end

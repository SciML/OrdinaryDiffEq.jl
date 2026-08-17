module OrdinaryDiffEqTaylorSeries

import OrdinaryDiffEqCore: alg_stability_size,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
    alg_cache,
    OrdinaryDiffEqConstantCache, trivial_limiter!,
    perform_step!, @cache,
    _ode_interpolant, _ode_interpolant!,
    OrdinaryDiffEqAlgorithm,
    _ode_addsteps!,
    get_fsalfirstlast, isfsal,
    DerivativeOrderNotPossibleError, unwrap_alg, step_accept_controller!,
    stepsize_controller!, get_current_adaptive_order, get_current_alg_order
using FastBroadcast: Serial
import MuladdMacro: @muladd
import FastBroadcast: @..
import RecursiveArrayTools: recursivefill!
using TruncatedStacktraces: @truncate_stacktrace
using TaylorDiff: TaylorDiff, TaylorArray, TaylorScalar
using Symbolics: Symbolics, @variables, build_function
using SymbolicUtils: SymbolicUtils
import CommonSolve: solve
import SciMLBase: SciMLBase, unwrapped_f, alg_order
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import OrdinaryDiffEqCore
import FunctionWrappers: FunctionWrapper

using Reexport: Reexport, @reexport
using SciMLBase: ODEFunction
@reexport using SciMLBase: ODEProblem, solve

include("algorithms.jl")
include("alg_utils.jl")
include("TaylorSeries_caches.jl")
include("TaylorSeries_perform_step.jl")
include("interpolants.jl")

import PrecompileTools
import Preferences

PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [ExplicitTaylor2()]
    prob_list = []

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

export ExplicitTaylor2, ExplicitTaylor, ExplicitTaylorAdaptiveOrder

end

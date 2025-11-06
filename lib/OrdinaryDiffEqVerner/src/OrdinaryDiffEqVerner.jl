module OrdinaryDiffEqVerner

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, unwrap_alg,
                           calculate_residuals, alg_stability_size,
                           OrdinaryDiffEqAlgorithm,
                           CompositeAlgorithm, AbstractController, PIDController,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val,
                           explicit_rk_docstring, trivial_limiter!, _ode_interpolant,
                           _ode_interpolant!, _ode_addsteps!, @fold,
                           @OnDemandTableauExtract, AutoAlgSwitch,
                           DerivativeOrderNotPossibleError,
                           get_fsalfirstlast, copyat_or_push!
using FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros
using Static: False
using TruncatedStacktraces: @truncate_stacktrace
using LinearAlgebra: norm
import OrdinaryDiffEqCore
using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("verner_tableaus.jl")
include("verner_caches.jl")
include("verner_addsteps.jl")
include("interp_func.jl")
include("interpolants.jl")
include("controllers.jl")
include("verner_rk_perform_step.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [Vern6(), Vern7(), Vern8(), Vern9()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]))
    end

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export Vern6, Vern7, Vern8, Vern9, RKV76IIa
export AutoVern6, AutoVern7, AutoVern8, AutoVern9

end

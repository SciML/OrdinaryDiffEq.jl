module OrdinaryDiffEqSSPRK

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, unwrap_alg,
                           calculate_residuals, ssp_coefficient,
                           OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                           OrdinaryDiffEqAdaptiveAlgorithm, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val,
                           explicit_rk_docstring, trivial_limiter!,
                           _ode_interpolant, _ode_interpolant!,
                           _ode_addsteps!, get_fsalfirstlast, @SciMLMessage, copyat_or_push!
using FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def
using Static: False
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

include("algorithms.jl")
include("alg_utils.jl")
include("ssprk_caches.jl")
include("interp_func.jl")
include("ssprk_perform_step.jl")
include("interpolants.jl")
include("addsteps.jl")
include("functions.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = []
    solver_list_nonadaptive = []
    prob_list = []

    low_storage = [
        SSPRK43()
    ]

    low_storage_nonadaptive = [
    ]

    if Preferences.@load_preference("PrecompileLowStorage", false)
        append!(solver_list, low_storage)
        append!(solver_list_nonadaptive, low_storage_nonadaptive)
    end

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

    for prob in prob_list, solver in solver_list_nonadaptive
        solve(prob, solver; dt = 0.5)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
       SSPRK33, KYKSSPRK42, KYK2014DGSSPRK_3S2

end

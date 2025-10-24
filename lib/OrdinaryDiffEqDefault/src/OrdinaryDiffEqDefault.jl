module OrdinaryDiffEqDefault

using OrdinaryDiffEqCore: alg_stability_size, beta2_default, beta1_default, AutoSwitchCache,
                          ODEIntegrator, trivial_limiter!,
                          CompositeAlgorithm, OrdinaryDiffEqAlgorithm,
                          OrdinaryDiffEqMutableCache, AutoAlgSwitch
using OrdinaryDiffEqVerner: Vern7, Vern8, Vern9, Vern6
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
using OrdinaryDiffEqBDF: FBDF
import OrdinaryDiffEqCore

import OrdinaryDiffEqCore: is_mass_matrix_alg, default_autoswitch, isdefaultalg, @SciMLMessage
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType
import LinearSolve
using LinearAlgebra: I, isdiag
using EnumX

using Reexport
@reexport using SciMLBase

include("default_alg.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = []
    prob_list = []

    default_ode = [
        DefaultODEAlgorithm(autodiff = AutoFiniteDiff())
    ]

    default_autodiff_ode = [
        DefaultODEAlgorithm()
    ]

    if Preferences.@load_preference("PrecompileDefault", true)
        append!(solver_list, default_ode)
    end

    if Preferences.@load_preference("PrecompileAutodiffDefault", true)
        append!(solver_list, default_autodiff_ode)
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

    prob_list = nothing
    solver_list = nothing
end

export DefaultODEAlgorithm, DefaultImplicitODEAlgorithm

end # module OrdinaryDiffEqDefault

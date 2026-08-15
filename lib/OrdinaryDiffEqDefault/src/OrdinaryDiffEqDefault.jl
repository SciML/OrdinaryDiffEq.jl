module OrdinaryDiffEqDefault

using OrdinaryDiffEqCore: alg_stability_size, beta2_default, beta1_default, AutoSwitchCache,
    CompositeAlgorithm, AutoAlgSwitch
using OrdinaryDiffEqVerner: Vern7
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
using OrdinaryDiffEqBDF: FBDF, DFBDF

import OrdinaryDiffEqCore: is_mass_matrix_alg, default_autoswitch, isdefaultalg
import ADTypes: AutoFiniteDiff
import DiffEqBase
import LinearSolve
using LinearAlgebra: I
using EnumX: EnumX

import SciMLBase

using Reexport: Reexport, @reexport
@reexport using SciMLBase: ODEProblem, ODEFunction, solve, init, solve!, step!, remake,
    reinit!, ReturnCode, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    CallbackSet, terminate!, DAEProblem, DAEFunction

include("default_alg.jl")

function _lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return
end

function _lorenz_p!(du, u, p, t)
    du[1] = p.σ * (u[2] - u[1])
    du[2] = u[1] * (p.ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - p.β * u[3]
    return
end

const _lorenz_p_params = (σ = 10.0, ρ = 28.0, β = 8 / 3)

function _lorenz_pref!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
    return
end

const _lorenz_pref_params = [10.0, 28.0, 8 / 3]

function _dae_p!(resid, du, u, p, t)
    resid[1] = du[1] + p.rate * u[1]
    return nothing
end

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    solver_list = []
    prob_list = []

    default_ode = [
        DefaultODEAlgorithm(autodiff = AutoFiniteDiff()),
    ]

    default_autodiff_ode = [
        DefaultODEAlgorithm(),
    ]

    if Preferences.@load_preference("PrecompileDefault", true)
        append!(solver_list, default_ode)
    end

    if Preferences.@load_preference("PrecompileAutodiffDefault", true)
        append!(solver_list, default_autodiff_ode)
    end

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(_lorenz!, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(_lorenz!, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    # Kept in their own list as well as `prob_list`: besides the usual
    # `solve(prob, solver)` sweep, these also drive the no-algorithm
    # `solve(prob)` entry below.
    parameter_generic_prob_list = []

    if Preferences.@load_preference("PrecompileAutoDespecialize", true)
        prob = ODEProblem{true, SciMLBase.AutoDespecialize}(
            _lorenz_p!, [1.0; 0.0; 0.0],
            (0.0, 1.0), _lorenz_p_params
        )
        push!(parameter_generic_prob_list, prob)

        dae_f = DAEFunction{true, SciMLBase.AutoDespecialize}(_dae_p!)
        dae_prob = DAEProblem(dae_f, [-0.5], [1.0], (0.0, 0.1), (rate = 0.5,))
        solve(dae_prob)
    end

    if Preferences.@load_preference("PrecompileAutoDePSpecialize", true)
        push!(
            parameter_generic_prob_list,
            ODEProblem{true, SciMLBase.AutoDePSpecialize}(
                _lorenz_p!, [1.0; 0.0; 0.0],
                (0.0, 1.0), _lorenz_p_params
            )
        )
        push!(
            parameter_generic_prob_list,
            ODEProblem{true, SciMLBase.AutoDePSpecialize}(
                _lorenz_pref!, [1.0; 0.0; 0.0],
                (0.0, 1.0), _lorenz_pref_params
            )
        )
    end

    append!(prob_list, parameter_generic_prob_list)

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0], (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(
                _lorenz!, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]
            )
        )
    end

    for prob in prob_list, solver in solver_list

        solve(prob, solver)(5.0)
    end

    # `solve(prob)` runs default algorithm selection, a separate entry point from
    # `solve(prob, DefaultODEAlgorithm())` above.
    for prob in parameter_generic_prob_list

        solve(prob)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
    parameter_generic_prob_list = nothing
end

export DefaultODEAlgorithm, DefaultImplicitODEAlgorithm

end # module OrdinaryDiffEqDefault

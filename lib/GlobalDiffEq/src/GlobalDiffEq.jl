module GlobalDiffEq

using Reexport: @reexport
@reexport using DiffEqBase

import ADTypes, DifferentiationInterface, ForwardDiff, LinearAlgebra,
    OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, Random, RecursiveArrayTools,
    Richardson, SciMLBase, SciMLStructures
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import OrdinaryDiffEqCore: perform_step!, @cache
using FastBroadcast: FastBroadcast, @..
using MuladdMacro: MuladdMacro, @muladd
using PrecompileTools: @setup_workload, @compile_workload

abstract type GlobalDiffEqAlgorithm <: SciMLBase.AbstractODEAlgorithm end

include("richardson.jl")
include("adjoint.jl")
include("companion.jl")
include("transport.jl")
include("glee/tableaus.jl")
include("glee/algorithms.jl")
include("glee/solve.jl")
include("glee/caches.jl")
include("glee/perform_step.jl")

export GlobalAdjoint, GlobalRichardson, adjoint_error_estimate
export GlobalErrorTransport
export GLEE23, GLEE24, GLEE35, MM5GEE, global_error_estimate

@setup_workload begin
    # Simple test ODE: exponential decay du/dt = -u
    function f!(du, u, p, t)
        du[1] = -u[1]
    end
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(f!, u0, tspan)

    @compile_workload begin
        solve(
            prob, GlobalRichardson(OrdinaryDiffEqTsit5.Tsit5()),
            dt = 0.1, reltol = 1.0e-3, abstol = 1.0e-6
        )
    end
end

end

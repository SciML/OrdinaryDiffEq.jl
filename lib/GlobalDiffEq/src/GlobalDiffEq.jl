module GlobalDiffEq

using Reexport: @reexport
@reexport using DiffEqBase

import OrdinaryDiffEqTsit5, Richardson, SciMLBase
using PrecompileTools: @setup_workload, @compile_workload

abstract type GlobalDiffEqAlgorithm <: SciMLBase.AbstractODEAlgorithm end

include("richardson.jl")

export GlobalRichardson

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

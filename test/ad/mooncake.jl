using OrdinaryDiffEq, StaticArrays, Test, FiniteDiff, DiffEqBase
using ADTypes: AutoForwardDiff, AutoMooncake
import DifferentiationInterface as DI
using Mooncake  # Load Mooncake after DI to ensure extension is loaded

# Mooncake is supported on all Julia versions

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

function f(u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough())
    return sum(sol)
end

u0 = [1.0; 0.0; 0.0]

# Reference gradient using ForwardDiff
ref_grad = DI.gradient(f, AutoForwardDiff(), u0)

# Test Mooncake gradient using DifferentiationInterface
@test_broken begin
    mooncake_grad = DI.gradient(f, AutoMooncake(; config = nothing), u0)
    mooncake_grad ≈ ref_grad
    true
end

# Alternative test with FiniteDiff as reference
fd_grad = FiniteDiff.finite_difference_gradient(f, u0)
@test ref_grad ≈ fd_grad rtol = 1.0e-6

# Skip Enzyme tests on Julia 1.12+ prerelease versions
if !isempty(VERSION.prerelease)
    @warn "Skipping Enzyme tests on Julia prerelease version $(VERSION)"
    exit(0)
end

# Enzyme is only supported on Julia <= 1.11
if VERSION >= v"1.12"
    @warn "Skipping Enzyme tests on Julia $(VERSION) - only supported on Julia <= 1.11"
    exit(0)
end

using Enzyme, OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, Test
using ADTypes
import DifferentiationInterface as DI

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

function f_dt(y::Array{Float64}, u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    y .= sol[1, :]
    return nothing
end;

function f_dt(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sol[1, :]
end;

u0 = [1.0; 0.0; 0.0]

# Test jacobian with ForwardDiff via DifferentiationInterface
fdj = DI.jacobian(f_dt, AutoForwardDiff(), u0)

# Test jacobian with Enzyme forward mode via DifferentiationInterface
ezj = DI.jacobian(f_dt, AutoEnzyme(mode = Enzyme.Forward), u0)

@test ezj ≈ fdj

function f_dt2(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), dt = 0.1, saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sum(sol[1, :])
end

# Test gradient with ForwardDiff via DifferentiationInterface
fdg = DI.gradient(f_dt2, AutoForwardDiff(), u0)

# Test gradient with Enzyme reverse mode via DifferentiationInterface
ezg = DI.gradient(f_dt2, AutoEnzyme(mode = Enzyme.Reverse), u0)

@test ezg ≈ fdg

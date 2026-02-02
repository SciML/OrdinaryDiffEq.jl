# Discrete adjoint tests with Mooncake
# Enzyme: skipped due to segfaults (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
# Mooncake: all versions (currently broken)
# ForwardDiff: all versions (reference)

using OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, Test, ForwardDiff
using ADTypes
import DifferentiationInterface as DI
using Mooncake  # Load Mooncake after DI to ensure extension is loaded

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

const _saveat = SA[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

# Functions for jacobian tests (output array)
function f_dt(y::Array{Float64}, u0::Array{Float64})
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    y .= sol[1, :]
    return nothing
end

function f_dt(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sol[1, :]
end

# Function for gradient tests (scalar output)
function f_dt_sum(u0)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz!, u0, tspan)
    sol = DiffEqBase.solve(prob, Tsit5(), dt = 0.1, saveat = _saveat, sensealg = DiffEqBase.SensitivityADPassThrough(), abstol = 1.0e-12, reltol = 1.0e-12)
    return sum(sol[1, :])
end

u0 = [1.0; 0.0; 0.0]

# Reference jacobian and gradient using ForwardDiff
fdj = DI.jacobian(f_dt, AutoForwardDiff(), u0)
fdg = DI.gradient(f_dt_sum, AutoForwardDiff(), u0)

@testset "Discrete Adjoints" begin
    # Enzyme tests skipped - Enzyme segfaults on ODE solves which crashes the process
    # before @test_broken can catch it. See https://github.com/EnzymeAD/Enzyme.jl/issues/2699

    # Mooncake tests (all Julia versions)
    @testset "Mooncake" begin
        @testset "Gradient (Reverse mode)" begin
            # Mooncake is a reverse-mode AD, so we test gradients
            @test_broken begin
                mkg = DI.gradient(f_dt_sum, AutoMooncake(; config = nothing), u0)
                mkg ≈ fdg
            end
        end
    end
end

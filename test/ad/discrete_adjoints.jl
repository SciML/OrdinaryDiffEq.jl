# Discrete adjoint tests with Enzyme and Mooncake
# Version-dependent AD backend selection
# Enzyme: Julia <= 1.11 only (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
# Mooncake: all versions
# ForwardDiff: all versions (reference)

using OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, Test
using ADTypes
import DifferentiationInterface as DI
using Mooncake  # Load Mooncake after DI to ensure extension is loaded

const JULIA_VERSION_ALLOWS_ENZYME = VERSION < v"1.12" && isempty(VERSION.prerelease)

if JULIA_VERSION_ALLOWS_ENZYME
    using Enzyme
end

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
    # Enzyme tests (Julia <= 1.11 only)
    if JULIA_VERSION_ALLOWS_ENZYME
        @testset "Enzyme" begin
            @testset "Jacobian (Forward mode)" begin
                ezj = DI.jacobian(f_dt, AutoEnzyme(mode = Enzyme.Forward), u0)
                @test ezj ≈ fdj
            end

            @testset "Gradient (Reverse mode)" begin
                ezg = DI.gradient(f_dt_sum, AutoEnzyme(mode = Enzyme.Reverse), u0)
                @test ezg ≈ fdg
            end
        end
    end

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

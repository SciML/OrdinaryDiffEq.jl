# Discrete adjoint tests with Mooncake
# Enzyme: skipped due to segfaults (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
# Mooncake: all versions
#   - SensitivityADPassThrough (true discrete adjoint): broken — Mooncake cannot
#     differentiate through the mutating solver internals and silently produces
#     incorrect gradients.
#   - Continuous adjoints via SciMLSensitivity (e.g. InterpolatingAdjoint with
#     ZygoteVJP/EnzymeVJP, GaussAdjoint): working on all Julia versions, tested
#     below as a regression check that the SciMLSensitivityMooncakeExt path works.
# ForwardDiff: all versions (reference)

using OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, Test, ForwardDiff
using SciMLSensitivity
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

# Continuous-adjoint variant (parameterised Lorenz, Vector saveat, no
# SensitivityADPassThrough) so the SciMLSensitivity adjoint dispatch is used.
# Mooncake supports this path via SciMLSensitivityMooncakeExt with
# InterpolatingAdjoint(ZygoteVJP/EnzymeVJP) and GaussAdjoint(ZygoteVJP).
# A parameterised form is required because the no-parameter `lorenz!` triggers
# a "nothing returned from Zygote vjp" error in ZygoteVJP-based sensealgs.
function lorenz_p!(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    return du[3] = u[1] * u[2] - p[3] * u[3]
end

const _saveat_v = collect(0.0:0.25:3.0)
const _lorenz_p = [10.0, 28.0, 8 / 3]

function f_dt_sum_sa(u0, sensealg)
    tspan = (0.0, 3.0)
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz_p!, u0, tspan, _lorenz_p)
    sol = solve(prob, Tsit5(), saveat = _saveat_v, sensealg = sensealg, abstol = 1.0e-12, reltol = 1.0e-12)
    return sum(sol[1, :])
end

u0 = [1.0; 0.0; 0.0]

# Reference jacobian and gradient using ForwardDiff
fdj = DI.jacobian(f_dt, AutoForwardDiff(), u0)
fdg = DI.gradient(f_dt_sum, AutoForwardDiff(), u0)
fdg_sa = DI.gradient(u_ -> f_dt_sum_sa(u_, ForwardDiffSensitivity()), AutoForwardDiff(), u0)

@testset "Discrete Adjoints" begin
    # Enzyme tests skipped - Enzyme segfaults on ODE solves which crashes the process
    # before @test_broken can catch it. See https://github.com/EnzymeAD/Enzyme.jl/issues/2699

    # Mooncake tests (all Julia versions)
    @testset "Mooncake" begin
        # SensitivityADPassThrough lets Mooncake try to differentiate through the
        # mutating solver internals directly. Mooncake currently does not handle
        # all the in-place mutations correctly and silently produces wrong
        # gradients (the call succeeds but the result is incorrect), so this is
        # left as @test_broken.
        @testset "Gradient via SensitivityADPassThrough (broken)" begin
            @test_broken begin
                mkg = DI.gradient(f_dt_sum, AutoMooncake(; config = nothing), u0)
                mkg ≈ fdg
            end
        end

        # Continuous adjoints through SciMLSensitivity DO work with Mooncake
        # via the SciMLSensitivityMooncakeExt extension. These cover the same
        # Lorenz problem as the SensitivityADPassThrough test above.
        @testset "Gradient via $(nameof(typeof(sensealg))) ($(nameof(typeof(sensealg.autojacvec))))" for sensealg in (
                InterpolatingAdjoint(autojacvec = ZygoteVJP()),
                InterpolatingAdjoint(autojacvec = EnzymeVJP()),
                GaussAdjoint(autojacvec = ZygoteVJP()),
            )
            mkg = DI.gradient(u_ -> f_dt_sum_sa(u_, sensealg), AutoMooncake(; config = nothing), u0)
            @test mkg ≈ fdg_sa rtol = 1.0e-6
        end
    end
end

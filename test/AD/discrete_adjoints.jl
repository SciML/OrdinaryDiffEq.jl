# Discrete adjoint tests with Mooncake
# Enzyme: skipped due to segfaults (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
# Mooncake: all versions
# ForwardDiff: all versions (reference)

using OrdinaryDiffEqTsit5, StaticArrays, DiffEqBase, Test, ForwardDiff
using ADTypes
import DifferentiationInterface as DI
using Mooncake  # Load Mooncake after DI to ensure extension is loaded

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return
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

f_dt_sum_scalar(x) = f_dt_sum([x, 0.0, 0.0])

u0 = [1.0; 0.0; 0.0]

# Reference jacobian and gradient using ForwardDiff
fdj = DI.jacobian(f_dt, AutoForwardDiff(), u0)
fdg = DI.gradient(f_dt_sum, AutoForwardDiff(), u0)
fdd = DI.derivative(f_dt_sum_scalar, AutoForwardDiff(), 1.0)
@testset "Discrete Adjoints" begin
    # Enzyme tests skipped - Enzyme segfaults on ODE solves which crashes the process
    # before @test_broken can catch it. See https://github.com/EnzymeAD/Enzyme.jl/issues/2699

    # Mooncake tests (all Julia versions)
    @testset "Mooncake" begin
        @testset "Gradient via SensitivityADPassThrough" begin
            mkg = DI.gradient(f_dt_sum, AutoMooncake(; config = nothing), u0)
            @test mkg isa typeof(fdg)
            @test mkg ≈ fdg rtol = 1.0e-6
        end
        @testset "Derivative via SensitivityADPassThrough" begin
            mkd = DI.derivative(f_dt_sum_scalar, AutoMooncake(; config = nothing), 1.0)
            @test mkd isa typeof(fdd)
            @test mkd ≈ fdd rtol = 1.0e-6
        end
        @testset "promote_f rdata over p types and specializations" begin
            # `rdata_type` on a tuple keeps the full 2-tuple unless *every* field is
            # `NoRData`, so the shape of the pullback's `dy` depends on `p`. The tests
            # above only exercise an array `p`, which is the one configuration where the
            # old unpacking happened to be right.
            lin!(du, u, p, t) = (du .= -0.3 .* u)
            linp!(du, u, p, t) = (du .= -(p isa Number ? p : p[1]) .* u)
            v0 = [1.0, 2.0]

            for p in (SciMLBase.NullParameters(), [0.3, 1.0], 0.3, (0.3, 0.5))
                for spec in (SciMLBase.FullSpecialize, SciMLBase.NoSpecialize)
                    rhs = p isa SciMLBase.NullParameters ? lin! : linp!
                    g = let p = p, spec = spec, rhs = rhs
                        function (x)
                            prob = ODEProblem{true, spec}(rhs, x, (0.0, 1.0), p)
                            sol = DiffEqBase.solve(
                                prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10,
                                sensealg = DiffEqBase.SensitivityADPassThrough()
                            )
                            return sum(abs2, sol.u[end])
                        end
                    end
                    reference = DI.gradient(g, AutoForwardDiff(), v0)
                    @test DI.gradient(g, AutoMooncake(; config = nothing), v0) ≈
                        reference rtol = 1.0e-6
                end
            end
        end
    end
end

using OrdinaryDiffEqBDF, StaticArrays, LinearAlgebra, Test, ADTypes
using OrdinaryDiffEqNonlinearSolve

f = function (du, u, p, t)
    out1 = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    out2 = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    out3 = u[1] + u[2] + u[3] - 1.0
    return [out1, out2, out3]
end

u₀ = [1.0, 0, 0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

integrator = init(prob, DImplicitEuler())

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

integrator = init(prob, DFBDF())

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
@test integrator.u ≈ [1.0, 0, 0.0] atol = 1.0e-9
integrator = init(
    prob, DABDF2(), initializealg = OrdinaryDiffEqNonlinearSolve.ShampineCollocationInit()
)
@test !(integrator.u ≈ [1.0, 0, 0.0])

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan)
integrator = init(prob, DABDF2())
@test !(integrator.u ≈ [1.0, 0, 0.0])

f = function (out, du, u, p, t)
    out[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    return out[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
integrator2 = init(prob, DABDF2(autodiff = AutoFiniteDiff()))

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

@test integrator2.du[1] ≈ -0.04 atol = 1.0e-99
@test integrator2.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator2.u ≈ u₀ atol = 1.0e-9

integrator = init(prob, DImplicitEuler())

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

integrator = init(prob, DFBDF())

@test integrator.du[1] ≈ -0.04 atol = 1.0e-9
@test integrator.du[2] ≈ 0.04 atol = 1.0e-9
@test integrator.u ≈ u₀ atol = 1.0e-9

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
@test integrator.u ≈ [1.0, 0, 0.0] atol = 1.0e-9
integrator = init(
    prob, DABDF2(), initializealg = OrdinaryDiffEqNonlinearSolve.ShampineCollocationInit()
)
@test !(integrator.u ≈ [1.0, 0, 0.0])

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan)
integrator = init(prob, DABDF2())
@test !(integrator.u ≈ [1.0, 0, 0.0])

# Need to be able to find the consistent solution of this problem, broken right now
# analytical solution:
#   u[1](t) ->  cos(t)
#   u[2](t) -> -sin(t)
#   u[3](t) -> 2cos(t)
f = function (out, du, u, p, t)
    out[1] = du[1] - u[2]
    out[2] = du[2] + u[3] - cos(t)
    return out[3] = u[1] - cos(t)
end

u₀ = [1.0, 0.0, 0.0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 1.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(
    prob, DABDF2(); initializealg = OrdinaryDiffEqNonlinearSolve.ShampineCollocationInit()
)

@test integrator.du[1] ≈ 0.0 atol = 1.0e-9
@test_broken integrator.du[2] ≈ -1.0 atol = 1.0e-9
@test_broken integrator.u[3] ≈ 2.0 atol = 1.0e-9

struct UnusedParam
end

# test iip dae initialization with parameters without eltype/length
probp = DAEProblem(f, du₀, u₀, tspan, UnusedParam(), differential_vars = differential_vars)
for initializealg in (
        OrdinaryDiffEqNonlinearSolve.ShampineCollocationInit(),
        OrdinaryDiffEqNonlinearSolve.BrownFullBasicInit(),
    )
    @test isapprox(
        init(probp, DABDF2(); initializealg).u, init(prob, DABDF2(); initializealg).u
    )
end

f = function (du, u, p, t)
    return du - u
end

u₀ = SVector(1.0)
du₀ = SVector(0.0)
tspan = (0.0, 1.0)
differential_vars = SVector(true)
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du ≈ [1.0] atol = 1.0e-9

f = function (du, u, p, t)
    return du .- u
end

u₀ = SA[1.0, 1.0]
du₀ = SA[0.0, 0.0]
tspan = (0.0, 1.0)
differential_vars = [true, true]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du[1] ≈ 1.0 atol = 1.0e-9
@test integrator.du[2] ≈ 1.0 atol = 1.0e-9
# test oop DAE initialization with parameters without eltype/length
probp = DAEProblem(f, du₀, u₀, tspan, UnusedParam(), differential_vars = differential_vars)
for initializealg in (
        OrdinaryDiffEqNonlinearSolve.ShampineCollocationInit(),
        OrdinaryDiffEqNonlinearSolve.BrownFullBasicInit(),
    )
    @test isapprox(
        init(probp, DABDF2(); initializealg).u, init(prob, DABDF2(); initializealg).u
    )
end

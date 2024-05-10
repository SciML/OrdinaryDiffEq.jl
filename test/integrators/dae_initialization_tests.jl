using OrdinaryDiffEq, StaticArrays, LinearAlgebra, Test

## Mass Matrix

function rober_oop(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    [du1, du2, du3]
end
M = [1.0 0 0
     0 1.0 0
     0 0 0]
f_oop = ODEFunction(rober_oop, mass_matrix = M)
prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rosenbrock23(autodiff = false), reltol = 1e-8, abstol = 1e-8)
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(prob_mm, Rosenbrock23(), reltol = 1e-8, abstol = 1e-8,
    initializealg = ShampineCollocationInit())
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rosenbrock23(), reltol = 1e-8, abstol = 1e-8)
@test sum(sol[1]) ≈ 1
@test sol[1] ≈ [1.0, 0.0, 0.0]
for alg in [Rosenbrock23(autodiff = false), Trapezoid()]
    local sol
    sol = solve(prob_mm, alg, reltol = 1e-8, abstol = 1e-8,
        initializealg = ShampineCollocationInit())
    @test sum(sol[1]) ≈ 1
end

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
M = [1.0 0 0
     0 1.0 0
     0 0 0]
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5(autodiff = false), reltol = 1e-8, abstol = 1e-8)
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8,
    initializealg = ShampineCollocationInit())
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

prob_mm = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8)
@test sum(sol[1]) ≈ 1
@test sol[1] ≈ [1.0, 0.0, 0.0]

for alg in [Rodas5(autodiff = false), Trapezoid()]
    local sol
    sol = solve(prob_mm, alg, reltol = 1e-8, abstol = 1e-8,
        initializealg = ShampineCollocationInit())
    @test sum(sol[1]) ≈ 1
end

function rober_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3e7, 1e4)
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end

function rober_oop_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3e7, 1e4)
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    [du1, du2, du3]
end

# test oop and iip ODE initialization with parameters without eltype/length
struct UnusedParam
end
for f in (
    ODEFunction(rober_no_p, mass_matrix = M), ODEFunction(rober_oop_no_p, mass_matrix = M))
    local prob, probp
    prob = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1e5))
    probp = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1e5), UnusedParam)
    for initializealg in (ShampineCollocationInit(), BrownFullBasicInit())
        isapprox(init(prob, Rodas5(), abstol = 1e-10; initializealg).u,
            init(prob, Rodas5(), abstol = 1e-10; initializealg).u)
    end
end

## DAEProblem

f = function (du, u, p, t)
    out1 = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    out2 = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out3 = u[1] + u[2] + u[3] - 1.0
    [out1, out2, out3]
end

u₀ = [1.0, 0, 0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

integrator = init(prob, DImplicitEuler())

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

integrator = init(prob, DFBDF())

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
@test integrator.u≈[1.0, 0, 0.0] atol=1e-9
integrator = init(prob, DABDF2(), initializealg = ShampineCollocationInit())
@test !(integrator.u ≈ [1.0, 0, 0.0])

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan)
integrator = init(prob, DABDF2())
@test !(integrator.u ≈ [1.0, 0, 0.0])

f = function (out, du, u, p, t)
    out[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
integrator2 = init(prob, DABDF2(autodiff = false))

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

@test integrator2.du[1]≈-0.04 atol=1e-99
@test integrator2.du[2]≈0.04 atol=1e-9
@test integrator2.u≈u₀ atol=1e-9

integrator = init(prob, DImplicitEuler())

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

integrator = init(prob, DFBDF())

@test integrator.du[1]≈-0.04 atol=1e-9
@test integrator.du[2]≈0.04 atol=1e-9
@test integrator.u≈u₀ atol=1e-9

u₀ = [1.0, 0, 0.2]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())
@test integrator.u≈[1.0, 0, 0.0] atol=1e-9
integrator = init(prob, DABDF2(), initializealg = ShampineCollocationInit())
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
    out[3] = u[1] - cos(t)
end

u₀ = [1.0, 0.0, 0.0]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 1.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2(); initializealg = ShampineCollocationInit())

@test integrator.du[1]≈0.0 atol=1e-9
@test_broken integrator.du[2]≈-1.0 atol=1e-9
@test_broken integrator.u[3]≈2.0 atol=1e-9

# test iip dae initialization with parameters without eltype/length
probp = DAEProblem(f, du₀, u₀, tspan, UnusedParam(), differential_vars = differential_vars)
for initializealg in (ShampineCollocationInit(), BrownFullBasicInit())
    @test isapprox(
        init(probp, DABDF2(); initializealg).u, init(prob, DABDF2(); initializealg).u)
end

f = function (du, u, p, t)
    du - u
end

u₀ = SVector(1.0)
du₀ = SVector(0.0)
tspan = (0.0, 1.0)
differential_vars = SVector(true)
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du≈[1.0] atol=1e-9

f = function (du, u, p, t)
    du .- u
end

u₀ = SA[1.0, 1.0]
du₀ = SA[0.0, 0.0]
tspan = (0.0, 1.0)
differential_vars = [true, true]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
integrator = init(prob, DABDF2())

@test integrator.du[1]≈1.0 atol=1e-9
@test integrator.du[2]≈1.0 atol=1e-9
# test oop DAE initialization with parameters without eltype/length
probp = DAEProblem(f, du₀, u₀, tspan, UnusedParam(), differential_vars = differential_vars)
for initializealg in (ShampineCollocationInit(), BrownFullBasicInit())
    @test isapprox(
        init(probp, DABDF2(); initializealg).u, init(prob, DABDF2(); initializealg).u)
end

# to test that we get the right NL solve we need a broken solver.
struct BrokenNLSolve <: SciMLBase.AbstractNonlinearAlgorithm
    BrokenNLSolve(; kwargs...) = new()
end
function SciMLBase.__solve(prob::NonlinearProblem,
        alg::BrokenNLSolve, args...;
        kwargs...)
    u = fill(reinterpret(Float64, 0xDEADBEEFDEADBEEF), 3)
    SciMLBase.build_solution(prob, alg, u, copy(u);
        retcode = ReturnCode.Success)
end
function f2(u, p, t)
    u
end
f = ODEFunction(f2, mass_matrix = Diagonal([1.0, 1.0, 0.0]))
prob = ODEProblem(f, ones(3), (0.0, 1.0))
integrator = init(prob, Rodas5P(),
    initializealg = ShampineCollocationInit(1.0, BrokenNLSolve()))
@test all(isequal(reinterpret(Float64, 0xDEADBEEFDEADBEEF)), integrator.u)

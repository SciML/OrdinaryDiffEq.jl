using OrdinaryDiffEq, StaticArrays, LinearAlgebra, Test, ADTypes

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
sol = solve(
    prob_mm, Rosenbrock23(autodiff = AutoFiniteDiff()), reltol = 1e-8, abstol = 1e-8)
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(prob_mm, Rosenbrock23(), reltol = 1e-8, abstol = 1e-8,
    initializealg = ShampineCollocationInit())
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rosenbrock23(), reltol = 1e-8, abstol = 1e-8)
@test sum(sol[1]) ≈ 1
@test sol[1] ≈ [1.0, 0.0, 0.0]
for alg in [Rosenbrock23(autodiff = AutoFiniteDiff()), Trapezoid()]
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
sol = solve(prob_mm, Rodas5(autodiff = AutoFiniteDiff()), reltol = 1e-8, abstol = 1e-8)
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8,
    initializealg = ShampineCollocationInit())
@test sol[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

prob_mm = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8)
@test sum(sol[1]) ≈ 1
@test sol[1] ≈ [1.0, 0.0, 0.0]

for alg in [Rodas5(autodiff = AutoFiniteDiff()), Trapezoid()]
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

@testset "OverrideInit for DAEProblem" begin
    function daerhs(du, u, p, t)
        return [u[1] * t + p, u[1]^2 - u[2]^2]
    end
    # unknowns are u[2], p, D(u[1]), D(u[2]). Parameters are u[1], t
    initprob = NonlinearProblem([1.0, 1.0, 1.0, 1.0], [1.0, 0.0]) do x, _p
        u2, p, du1, du2 = x
        u1, t = _p
        return [u1^3 - u2^3, p^2 - 2p + 1, du1 - u1 * t - p, 2u1 * du1 - 2u2 * du2]
    end

    update_initializeprob! = function (iprob, integ)
        iprob.p[1] = integ.u[1]
        iprob.p[2] = integ.t
    end
    initprobmap = function (nlsol)
        return [parameter_values(nlsol)[1], nlsol.u[1]]
    end
    initprobpmap = function (_, nlsol)
        return nlsol.u[2]
    end
    initprob_du0map = function (nlsol)
        return nlsol.u[3:4]
    end
    initialization_data = SciMLBase.OverrideInitData(
        initprob, update_initializeprob!, initprobmap, initprobpmap, initprob_du0map)
    fn = DAEFunction(daerhs; initialization_data)
    prob = DAEProblem(fn, [0.0, 0.0], [2.0, 0.0], (0.0, 1.0), 0.0)
    integ = init(prob, DImplicitEuler())
    @test integ.du ≈ [1.0, 1.0]
    @test integ.u ≈ [2.0, 2.0]
    @test integ.p ≈ 1.0
    @test integ.sol.retcode != SciMLBase.ReturnCode.InitialFailure
end

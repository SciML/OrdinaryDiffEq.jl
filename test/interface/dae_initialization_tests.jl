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

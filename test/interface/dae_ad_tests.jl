using OrdinaryDiffEq, ForwardDiff, Test

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    [-k₁ * y₁ + k₃ * y₂ * y₃,
      k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
      y₁ + y₂ + y₃ - 1]
end
M = [1.0 0 0
    0 1.0 0
    0 0 0]
roberf = ODEFunction(rober, mass_matrix = M)
roberf_oop = ODEFunction{false}(rober, mass_matrix = M)
prob_mm = ODEProblem(roberf, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
prob_mm_oop = ODEProblem(roberf_oop, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5P(), reltol = 1e-8, abstol = 1e-8)
sol = solve(prob_mm_oop, Rodas5P(), reltol = 1e-8, abstol = 1e-8)

function f(out, du, u, p, t)
    out[1] = -p[1] * u[1] + p[3] * u[2] * u[3] - du[1]
    out[2] = +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end
function f(du, u, p, t)
    [-p[1] * u[1] + p[3] * u[2] * u[3] - du[1],
     +p[1] * u[1] - p[2] * u[2]^2 - p[3] * u[2] * u[3] - du[2],
     u[1] + u[2] + u[3] - 1.0]
end
p = [0.04, 3e7, 1e4]
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
prob_oop = DAEProblem{false}(f, du₀, u₀, tspan, p, differential_vars = differential_vars)
sol1 = solve(prob, DFBDF(), dt=1e-5, abstol = 1e-8, reltol = 1e-8)

# These tests flex differentiation of the solver and through the initialization
# To only test the solver part and isolate potential issues, set the initialization to consistent
@testset "Inplace: $(isinplace(_prob)), DAEProblem: $(_prob isa DAEProblem), BrownBasic: $(initalg isa BrownFullBasicInit)" for _prob in [prob, prob_mm, prob_oop, prob_mm_oop], initalg in [BrownFullBasicInit(), ShampineCollocationInit()]
    alg = _prob isa DAEProblem ? DFBDF() : Rodas5P()
    function f(p)
        sol = solve(remake(_prob, p = p), alg, abstol = 1e-14, reltol = 1e-14, initializealg = initalg)
        sum(sol)
    end
    @test ForwardDiff.gradient(f, [0.04, 3e7, 1e4]) ≈ [0,0,0] atol=1e-8
end

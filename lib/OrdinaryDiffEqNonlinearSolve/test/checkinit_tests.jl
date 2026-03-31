using OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock, LinearAlgebra, ForwardDiff, Test
using OrdinaryDiffEqCore

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    return [
        -k₁ * y₁ + k₃ * y₂ * y₃,
        k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
        y₁ + y₂ + y₃ - 1,
    ]
end
M = [
    1.0 0 0
    0 1.0 0
    0 0 0
]
roberf = ODEFunction(rober, mass_matrix = M)
roberf_oop = ODEFunction{false}(rober, mass_matrix = M)
prob_mm = ODEProblem(roberf, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
prob_mm_oop = ODEProblem(roberf_oop, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))

@test_throws SciMLBase.CheckInitFailureError solve(
    prob_mm, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8, initializealg = SciMLBase.CheckInit()
)
@test_throws SciMLBase.CheckInitFailureError solve(
    prob_mm_oop, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = SciMLBase.CheckInit()
)

f_oop = function (du, u, p, t)
    out1 = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    out2 = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    out3 = u[1] + u[2] + u[3] - 1.0
    return [out1, out2, out3]
end

f = function (resid, du, u, p, t)
    resid[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
    resid[2] = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
    return resid[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0.2]
du₀ = [0.0, 0.0, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)
prob_oop = DAEProblem(f_oop, du₀, u₀, tspan, differential_vars = differential_vars)
@test_throws SciMLBase.CheckInitFailureError solve(
    prob, DFBDF(), reltol = 1.0e-8, abstol = 1.0e-8, initializealg = SciMLBase.CheckInit()
)
@test_throws SciMLBase.CheckInitFailureError solve(
    prob_oop, DFBDF(), reltol = 1.0e-8, abstol = 1.0e-8, initializealg = SciMLBase.CheckInit()
)

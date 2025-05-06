using OrdinaryDiffEq, DiffEqDevTools, Test
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

prob = prob_ode_linear
sol = solve(prob, Rosenbrock32())
dt₀ = sol.t[2]

prob = prob_ode_2Dlinear
sol = solve(prob, ExplicitRK(tableau = constructBogakiShampine3()))
dt₀ = sol.t[2]

@test 1e-7 < dt₀ < 0.1
@test_throws ArgumentError local sol = solve(prob, Euler())
#dt₀ = sol.t[2]

sol3 = solve(prob, ExplicitRK(tableau = constructDormandPrince8_64bit()))
dt₀ = sol3.t[2]

@test 1e-7 < dt₀ < 0.3

T = Float32
u0 = T.([1.0; 0.0; 0.0])

tspan = T.((0, 70))
prob = remake(prob, u0 = u0, tspan = tspan)
@test_nowarn solve(prob, Euler(); dt = T(0.0001))

tspan = T.((2000, 2100))
prob = remake(prob, tspan = tspan)
# set maxiters to prevent infinite loop on test failure
@test solve(prob, Euler(); dt = T(0.0001), maxiters = 10).retcode ==
      SciMLBase.ReturnCode.MaxIters

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    nothing
end
u0 = Float32[1.0, 0.0, 0.0]
tspan = (0.0f0, 1.0f5)
params = (4.0f-2, 3.0f7, 1.0f4)
prob = ODEProblem(rober, u0, tspan, params)
sol = solve(prob, Rosenbrock23())

# https://github.com/SciML/DifferentialEquations.jl/issues/743

using LinearAlgebra
function f(du, u, p, t)
    du[1] = -p[1] * u[1] + p[2] * u[2] * u[3]
    du[2] = p[1] * u[1] - p[2] * u[2] * u[3] - p[3] * u[2] * u[2]
    du[3] = u[1] + u[2] + u[3] - 1.0
end
M = Diagonal([1, 1, 0])
p = [0.04, 10^4, 3e7]
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 1e6)
prob = ODEProblem(ODEFunction(f, mass_matrix = M), u0, tspan, p)
sol = solve(prob, Rodas5())
@test sol.t[end] == 1e6

# test that dtmin is set based on timespan
prob = ODEProblem((u, p, t) -> 1e20 * sin(1e20 * t), 0.1, (0, 1e-19))
@test solve(prob, Tsit5()).retcode == ReturnCode.Success

#test that we are robust to u0=0, t0!=0
integ = init(ODEProblem(((u, p, t) -> u), 0.0f0, (20.0f0, 0.0f0)), Tsit5())
@test abs(integ.dt) > eps(integ.t)
integ = init(ODEProblem(((du, u, p, t) -> du .= u), [0.0f0], (20.0f0, 0.0f0)), Tsit5())
@test abs(integ.dt) > eps(integ.t)

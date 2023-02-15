using OrdinaryDiffEq, DiffEqDevTools, Test

import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

prob = prob_ode_2Dlinear
sol = solve(prob, Rosenbrock32(), dt = 1 / 2^4)

sol = solve(prob, ExplicitRK(tableau = constructBogakiShampine3()), dt = 1 / 2^4)
val1 = maximum(abs.(sol.u[end] - sol.u_analytic[end]))

sol2 = solve(prob, ExplicitRK(tableau = constructDormandPrince()), dt = 1 / 2^4)
val2 = maximum(abs.(sol2.u[end] - sol2.u_analytic[end]))

sol3 = solve(prob, ExplicitRK(tableau = constructRKF8(Float64)), dt = 1 / 2^4)
val3 = maximum(abs.(sol3.u[end] - sol3.u_analytic[end]))

@test length(sol.t) > length(sol2.t) >= length(sol3.t)
@test max(val1, val2, val3) < 2e-3

function lorenz(u, p, t)
    [10.0(u[2] - u[1])
     u[1] * (28.0 - u[3]) - u[2]
     u[1] * u[2] - (8 / 3) * u[3]]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem{false}(lorenz, u0, tspan)
sol = solve(prob, QNDF())
@test length(sol.t) < 5000
sol = solve(prob, FBDF())
@test length(sol.t) < 6600

function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem{true}(lorenz, u0, tspan)
sol = solve(prob, QNDF())
@test length(sol.t) < 5000
sol = solve(prob, FBDF())
@test length(sol.t) < 6600

function lorenz(out, du, u, p, t)
    out[1] = 10.0(u[2] - u[1]) - du[1]
    out[2] = u[1] * (28.0 - u[3]) - u[2] - du[2]
    out[3] = u[1] * u[2] - (8 / 3) * u[3] - du[3]
end
u0 = [1.0; 0.0; 0.0]
du0 = [0.0; 0.0; 0.0]
tspan = (0.0, 100.0)
differential_vars = [true, true, true]
prob = DAEProblem(lorenz, du0, u0, tspan, differential_vars = differential_vars)
sol = solve(prob, DFBDF())
@test length(sol.t) < 6600

function lorenz(du, u, p, t)
    [10.0(u[2] - u[1]) - du[1]
     u[1] * (28.0 - u[3]) - u[2] - du[2]
     u[1] * u[2] - (8 / 3) * u[3] - du[3]]
end
u0 = [1.0; 0.0; 0.0]
du0 = [0.0; 0.0; 0.0]
tspan = (0.0, 100.0)
differential_vars = [true, true, true]
prob = DAEProblem{false}(lorenz, du0, u0, tspan, differential_vars = differential_vars)
sol = solve(prob, DFBDF())
@test length(sol.t) < 6600

rr(x1, x2) = (x1 * (-2.1474936f0 * (x2 + x1)))
possibly_singular(u, p, t) = [-rr(u...), rr(u...)]
tspan = (1.6078221f0, 2.0f0)
initial_condition = [2.1349438f6, -2.1349438f6]
prob = ODEProblem(possibly_singular, initial_condition, tspan)
for alg in [
    Rosenbrock23(),
    Rosenbrock32(),
    Rodas3(),
    Rodas4(),
    Rodas4P(),
    Rodas5(),
    Rodas5P(),
]
    sol = solve(prob, alg)
    @test sol.retcode == ReturnCode.Success
end

# test problems with zero-length vectors
ode = ODEProblem((du, u, semi, t) -> du .= u, Float64[], (0.0, 1.0))
@test_nowarn solve(ode, Tsit5())

# test if the early exit in ode_determine_initdt doesn't hit in the case of
# zero-length vectors, see https://github.com/SciML/OrdinaryDiffEq.jl/pull/1865
integrator = init(ode, Tsit5())
@test integrator.dt ≈ 1.0e-6

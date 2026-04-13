using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test,
    ODEInterface, ODEInterfaceDiffEq

import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_linear,
    prob_ode_2Dlinear,
    prob_ode_bigfloat2Dlinear

probbig = prob_ode_bigfloat2Dlinear
probnum = prob_ode_linear
probnumbig = prob_ode_bigfloatlinear
prob = prob_ode_2Dlinear

dts = (1 / 2) .^ (7:-1:4)
testTol = 0.2
bools = Vector{Bool}(undef, 0)

## DP5()
sim = test_convergence(dts, probnum, DP5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol
sim = test_convergence(dts, prob, DP5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol

tabalg = ExplicitRK()
sol1 = solve(probnum, DP5(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnum, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1.0e-10

sol1 = solve(prob, DP5(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(prob, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 3.0e-10)

sol1 = solve(probnum, DP5(), dt = 1 / 2^6, controller = PIController(0.14, 0.04))
sol2 = solve(probnum, tabalg, dt = 1 / 2^6, controller = PIController(0.14, 0.04))

# Should be identical
sol1 = solve(probnum, DP5())
sol2 = solve(probnum, tabalg, controller = PIController(0.17, 0.04))
sol3 = solve(probnum, dopri5())

@test sol1.t ≈ sol2.t
@test sol1.t ≈ sol3.t atol = 1.0e-6

sol1 = solve(prob, DP5(), dt = 1 / 8)
sol2 = solve(prob, tabalg, controller = PIController(0.17, 0.04), dt = 1 / 8)
sol3 = solve(prob, dopri5(), dt = 1 / 8)

@test sol1.t ≈ sol2.t
@test sol1.t ≈ sol3.t atol = 1.0e-5

sol4 = solve(prob, DP5(), dt = 1 / 8, calck = false)

@test sol1.t == sol4.t

### DP8()

dts = (1 / 2) .^ (3:-1:1)
sim = test_convergence(dts, probnumbig, DP8())
@test abs.(sim.𝒪est[:l2] - 8) < testTol
sim = test_convergence(dts, probbig, DP8())
@test abs.(sim.𝒪est[:l2] - 8) < testTol

sol1 = solve(probnum, DP8(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnum, DP8(), dt = 1 / 2^6)

# DP8() uses qmax_first_step=10000 on the first step (Sundials CVODE behavior)
# while dop853() (Fortran) does not, so adaptive solves may differ slightly
sol1 = solve(probnum, DP8())
sol2 = solve(probnum, dop853())

@test sol1.u[end] ≈ sol2.u[end] atol = 1.0e-6

sol1 = solve(prob, DP8(), dt = 1 / 2^6)
sol2 = solve(prob, dop853(), dt = 1 / 2^6)

@test sol1.u[end] ≈ sol2.u[end] atol = 1.0e-6

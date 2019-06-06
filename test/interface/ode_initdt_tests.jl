using OrdinaryDiffEq, DiffEqDevTools, Test
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

prob = prob_ode_linear
sol =solve(prob,Rosenbrock32())
dt₀ = sol.t[2]

prob = prob_ode_2Dlinear
sol =solve(prob,ExplicitRK(tableau=constructBogakiShampine3()))
dt₀ = sol.t[2]

@test  1e-7 < dt₀ < .1
@test_throws ErrorException sol = solve(prob,Euler())
#dt₀ = sol.t[2]

sol3 =solve(prob,ExplicitRK(tableau=constructDormandPrince8_64bit()))
dt₀ = sol3.t[2]

@test 1e-7 < dt₀ < .3

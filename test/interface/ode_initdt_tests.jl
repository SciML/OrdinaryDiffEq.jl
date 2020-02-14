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

T = Float32
u0 = T.([1.0;0.0;0.0])

tspan = T.((0,100))
prob = remake(prob, u0=u0, tspan=tspan)
@test_nowarn solve(prob, Euler(); dt=T(0.0001))

tspan = T.((2000,2100))
prob = remake(prob, tspan=tspan)
@test_throws ArgumentError solve(prob, Euler(); dt=T(0.0001)) # Loops forever

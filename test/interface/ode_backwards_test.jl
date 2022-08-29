using OrdinaryDiffEq, Test, Statistics
import ODEProblemLibrary: prob_ode_2Dlinear
prob = deepcopy(prob_ode_2Dlinear)
prob2 = ODEProblem(prob.f, prob.u0, (1.0, 0.0), 1.01)

sol = solve(prob2, DP5(), dt = -1 / 4, tstops = [0.5])
@test sol.t == [1.0, 0.75, 0.5, 0]

sol = solve(prob2, RK4(), dt = -1 / 4, adaptive = false)
@test sol.t == [1.00, 0.75, 0.5, 0.25, 0]

sol = solve(prob2, RK4(), dt = -1 / 4, tstops = [0.5, 0.33], adaptive = false)
@test sol.t ≈ [1.00, 0.75, 0.5, 0.33, 0.08, 0]

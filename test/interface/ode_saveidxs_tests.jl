using OrdinaryDiffEq, Test
import ODEProblemLibrary: prob_ode_linear,
    linear, f_2dlinear

# scalar, not in-place
prob = prob_ode_linear
sol = solve(prob, DP5(); save_idxs = 1)
sol(0.5) # test interpolation

# vector, not in-place
prob2 = ODEProblem(linear, [0.5], (0.0, 1.0), 1.01)
sol2 = solve(prob2, DP5(); save_idxs = 1)
@test sol.t ≈ sol2.t && sol.u ≈ sol2.u
sol2(0.5)

sol2b = solve(prob2, DP5(); save_idxs = [1])
@test sol.t ≈ sol2b.t && sol.u ≈ [u[1] for u in sol2b.u]
sol2b(0.5)

# vector, in-place
prob3 = ODEProblem(f_2dlinear, [0.5], (0.0, 1.0), 1.01)
sol3 = solve(prob3, DP5(); save_idxs = 1)
@test sol.t ≈ sol3.t && sol.u ≈ sol3.u
sol3(0.5)

sol3b = solve(prob3, DP5(); save_idxs = [1])
@test sol.t ≈ sol3b.t && sol2b.u ≈ sol3b.u
sol3b(0.5)

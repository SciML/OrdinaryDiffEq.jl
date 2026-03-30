using OrdinaryDiffEq, DiffEqBase, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

f = (u, p, t) -> u
fun = ODEFunction(f; analytic = (u0, p, t) -> u0 * exp(t))
prob = ODEProblem(fun, 1 / 2, (0.0, 1.0))

sol = solve(prob, Euler(); dt = 1 // 2^(4), dense_errors = true)
sol2 = solve(prob, Vern9(); dt = 1 // 2^(10), abstol = 1e-14, reltol = 1e-14)

prob2 = ODEProblem(f, 1 / 2, (0.0, 1.0))
sol3 = solve(prob_ode_linear, Euler(); dt = 1 // 2^(4))

errsol1 = appxtrue(sol, sol2)

sol4 = solve(prob, Euler(); dt = 1 // 2^(4))
test_sol = TestSolution(sol2)
errsol2 = appxtrue(sol4, test_sol)

sol5 = solve(prob, Euler(); dt = 1 // 2^(4))
test_sol = TestSolution(sol2.t, sol2.u[end])
errsol3 = appxtrue(sol5, test_sol)

@test errsol1.errors[:L2] ≈ 0.018865798306718855
@test errsol1.errors[:L2] ≈ errsol2.errors[:L2]

@test errsol1.errors[:final] ≈ sol.errors[:final]
@test errsol1.errors[:l2] ≈ sol.errors[:l2]
@test errsol1.errors[:l∞] ≈ sol.errors[:l∞]
@test errsol1.errors[:L∞] ≈ sol.errors[:L∞]
@test errsol1.errors[:L2] ≈ sol.errors[:L2]

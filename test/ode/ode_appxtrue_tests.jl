using OrdinaryDiffEq, DiffEqBase, DiffEqProblemLibrary, DiffEqDevTools

f = (t,u) -> u
prob = ODEProblem(f,1/2,[0;1.0])
analytic = (t,u0) -> u0*exp(t)

sol =solve(prob,Euler;dt=1//2^(4))

sol2 =solve(prob,Vern9;dt=1//2^(10))

prob2 = ODETestProblem(f,1/2,analytic,[0;1.0])
sol3 =solve(prob_ode_linear,Euler;dt=1//2^(4))

appxtrue!(sol,sol2)

sol4 =solve(prob,Euler;dt=1//2^(4))
test_sol = TestSolution(sol2)
appxtrue!(sol4,test_sol)

sol5 =solve(prob,Euler;dt=1//2^(4))
test_sol = TestSolution(sol2[end])
appxtrue!(sol5,test_sol)

sol.errors[:L2] ≈ 0.018865798306718855 && sol.errors[:L2] ≈ sol4.errors[:L2]

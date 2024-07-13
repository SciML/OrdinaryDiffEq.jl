using OrdinaryDiffEq, DiffEqDevTools,  Test, BenchmarkTools, LinearAlgebra


f = (u, p, t) -> [-u[2],u[1]]
prob = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [cos(t), sin(t)]),
                [1.0, 0.0],
                (0.0, 1.0))
invariant(x) = norm(x)

r = Relaxation(invariant)
sol_relax = solve(prob, Tsit5(); relaxation = r)

sol = solve(prob, Tsit5())
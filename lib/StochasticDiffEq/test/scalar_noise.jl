using StochasticDiffEq, DiffEqNoiseProcess, Test

f(du, u, p, t) = (du .= u)
g(du, u, p, t) = (du .= u)
u0 = rand(4, 2)

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRIW1())
sol = solve(prob, SRIW1())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRI())
sol = solve(prob, SRI())

W = WienerProcess(0.0, 0.0, 0.0, reset = false)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRI())
@test_throws ErrorException sol = solve(prob, SRI())

g(du, u, p, t) = (du .= 1)

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRA1())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRA2())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRA3())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SOSRA())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SOSRA2())

W = WienerProcess(0.0, 0.0, 0.0)
prob = SDEProblem(f, g, u0, (0.0, 1.0), noise = W)
sol = solve(prob, SRA())

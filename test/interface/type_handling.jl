using OrdinaryDiffEq
prob = ODEProblem((u, p, t) -> -u, BigFloat(1.0), (0.0, 1.0))
solve(prob, Tsit5())
solve(prob, KenCarp4())

# Function initial condition
f(u, p, t) = u
u0(p, t0) = ones(2)
prob = ODEProblem(f, u0, (0.0, 1.0))
sol = solve(prob, Tsit5())

# test different u and t type
ode_f(du, u, p, t) = du[1] = -u[1]
# autospecialize broken for the next 35 minutes
prob = ODEProblem(ODEFunction{true, SciMLBase.FullSpecialize()}(ode_f), [1f0], (0.0,1.0))
for alg in [Tsit5(), Rodas5P(), FBDF()]
    solve(prob, alg)
end

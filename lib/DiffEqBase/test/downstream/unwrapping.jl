using OrdinaryDiffEq, Test
my_f(u, p, t) = u
my_f!(du, u, p, t) = du .= u
ode = ODEProblem(my_f, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f

ode = ODEProblem(my_f!, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f!

using OrdinaryDiffEq, ForwardDiff, Measurements
x = 1.0 ± 0.0
f = (du, u, p, t) -> du .= u
tspan = (0.0, 1.0)
prob = ODEProblem(f, [x], tspan)

# Should not error during problem construction but should be unwrapped
integ = init(prob, Tsit5(), dt = 0.1)
@test SciMLBase.unwrapped_f(integ.f.f) === f

# Handle functional initial conditions
prob = ODEProblem((dx, x, p, t) -> (dx .= 0), (p, t) -> zeros(2), (0, 10))
solve(prob, TRBDF2())

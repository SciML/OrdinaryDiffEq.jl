using OrdinaryDiffEq
srand(100)

linear = (t,u) -> (1.01*u)
analytic_linear = (t,u0) -> u0*exp(1.01*t)
prob = ODEProblem(linear,1/2,[0,1//2,1],analytic=analytic_linear)

sol3 =solve(prob::ODEProblem,dt=1//2^(6))

1//2 ∈ sol3.t

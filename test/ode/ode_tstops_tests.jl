using OrdinaryDiffEq, DiffEqBase
srand(100)

linear = (t,u) -> (1.01*u)
analytic_linear = (t,u0) -> u0*exp(1.01*t)
prob = ODETestProblem(linear,1/2,analytic_linear,(0.0,1.0))

sol3 =solve(prob,Tsit5(),dt=1//2^(6),tstops=[1/2])

1//2 ∈ sol3.t

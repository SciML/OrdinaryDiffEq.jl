using OrdinaryDiffEq, Plots, DiffEqBase

## Start on Number
f = (t,u) -> (2u)
analytic = (t,u0) -> u0*exp(t)
prob = ODETestProblem(f,1/2+(1/4)im,analytic)

sol = solve(prob,Tsit5(),dt=1/2^4)

u0 = rand(Complex64,5,5,5)
prob = ODETestProblem(f,u0,analytic)

sol = solve(prob,Tsit5(),dt=1/2^4)

typeof(sol[1]) == Array{Complex64,3}

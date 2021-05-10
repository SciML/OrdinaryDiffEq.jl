using OrdinaryDiffEq
f(u,p,t) = 10.1*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, QNDF())
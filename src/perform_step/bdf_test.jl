using OrdinaryDiffEq
f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, QNDF())

f(u,p,t) = -10*t
u0= 1.0
tspan = (0.0,2.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, QNDF())

f(u,p,t) = 1.01*(u^0.7+u^(0.2)-sin(u))
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, QNDF())
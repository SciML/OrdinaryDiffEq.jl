using OrdinaryDiffEq

function lorenz(u,p,t)
 [10.0(u[2]-u[1])
  u[1]*(28.0-u[3]) - u[2]
 u[1]*u[2] - (8/3)*u[3]]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem{false}(lorenz,u0,tspan)
sol = solve(prob,QNDF(),dt=1e-6)



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
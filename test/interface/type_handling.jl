using OrdinaryDiffEq
prob = ODEProblem((u,p,t) -> -u,BigFloat(1.0),(0.0,1.0))
solve(prob,Tsit5())
solve(prob,KenCarp4())

# Function initial condition
f(u,p,t) = u
u0(p,t0) = ones(2)
prob = ODEProblem(f,u0,(0.0,1.0))
sol = solve(prob,Tsit5())

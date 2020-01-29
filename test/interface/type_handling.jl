using OrdinaryDiffEq
prob = ODEProblem((u,p,t) -> -u,BigFloat(1.0),(0.0,1.0))
solve(prob,Tsit5())
solve(prob,KenCarp4())

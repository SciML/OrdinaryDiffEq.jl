using OrdinaryDiffEq
u0=rand(300,20).*ones(300,20)/2
prob = prob_ode_2Dlinear_notinplace
prob2 = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,Euler(),dt=1//2^(4),save_timeseries=false)
sol =solve(prob2::ODEProblem,Euler(),dt=1//2^(4),save_timeseries=false)

@time sol =solve(prob::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)
@time sol2 =solve(prob2::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)

alloc1 = @allocated sol =solve(prob::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)
alloc2 = @allocated sol2 =solve(prob2::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)

alloc1 = @allocated sol =solve(prob::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)
alloc2 = @allocated sol2 =solve(prob2::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=false)

bool1 = alloc2 <= alloc1

sol = solve(prob_ode_large2Dlinear::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=true)
sol2 = solve(prob_ode_large2Dlinear::ODEProblem,Euler(),sol.timeseries,sol.t,sol.k;dt=1//2^(8),save_timeseries=true)

sol = solve(prob_ode_large2Dlinear::ODEProblem,Euler(),dt=1//2^(6),save_timeseries=true)
alloc1 = @allocated sol = solve(prob_ode_large2Dlinear::ODEProblem,Euler(),dt=1//2^(8),save_timeseries=true)
alloc2 = @allocated sol2 = solve(prob_ode_large2Dlinear::ODEProblem,Euler(),sol.timeseries,sol.t,sol.k;dt=1//2^(8),save_timeseries=true)

bool2 = alloc2 <= alloc1

bool1 && bool2

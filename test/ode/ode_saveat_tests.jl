using OrdinaryDiffEq, Plots

prob = prob_ode_linear

sol =solve(prob,DP5(),dt=1//2^(2),save_timeseries=false,dense=false)
sol2=solve(prob,DP5(),dt=1//2^(2),save_timeseries=false,dense=false,saveat=[1/2])

@test symdiff(sol.t,sol2.t) == [1/2]

sol =solve(prob,DP5(),dt=1//2^(2),save_timeseries=true,dense=true)
sol2=solve(prob,DP5(),dt=1//2^(2),save_timeseries=true,dense=true,saveat=[1/2])

sol2(.49)

interpd = sol2(0:1//2^(4):1)

#plot(0:1//2^(4):1,interpd)

@test symdiff(sol.t,sol2.t) == [1/2]

sol =solve(prob,RK4(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,RK4(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

sol =solve(prob,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

sol =solve(prob,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

prob = prob_ode_2Dlinear

sol =solve(prob,DP5(),dt=1//2^(2),save_timeseries=false,dense=true)
sol2=solve(prob,DP5(),dt=1//2^(2),save_timeseries=false,dense=true,saveat=[1/2])

@test symdiff(sol.t,sol2.t) == [1/2]

sol =solve(prob,DP5(),dt=1//2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,DP5(),dt=1//2^(2),save_timeseries=true,dense=false,saveat=[1/2])

@test symdiff(sol.t,sol2.t) == [1/2]

sol =solve(prob,RK4(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,RK4(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

sol =solve(prob,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

sol =solve(prob,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

@test symdiff(sol.t,sol2.t) == [.125,.6,.61,.8]

sol=solve(prob,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[0,.125,.6,.61,.8])

@test !(sol.t[2] ≈ 0)

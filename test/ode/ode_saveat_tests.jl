using OrdinaryDiffEq, Plots

bools = Vector{Bool}(0)
prob = prob_ode_linear

sol =solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=false,dense=false)
sol2=solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=false,dense=false,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=true,dense=true)
sol2=solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=true,dense=true,saveat=[1/2])

sol2(.49)

interpd = sol2(0:1//2^(4):1)

#plot(0:1//2^(4):1,interpd)

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,RK4(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,RK4(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=false,dense=true)
sol2=solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=false,dense=true,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,DP5(),dt=1//2^(2),save_timeseries=true,dense=false,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,RK4(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,RK4(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,Rosenbrock32(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false)
sol2=solve(prob::ODEProblem,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol=solve(prob::ODEProblem,Trapezoid(),dt=1/2^(2),save_timeseries=true,dense=false,saveat=[0,.125,.6,.61,.8])

push!(bools,!(sol.t[2] ≈ 0))

minimum(bools)

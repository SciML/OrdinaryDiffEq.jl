using OrdinaryDiffEq, RecursiveArrayTools, DiffEqBase, Base.Test#, ParameterizedFunctions


f = function (t,u)
  - u + sin(t)
end


prob = ODEProblem(f,1.0,(0.0,10.0))

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  t - 2.95
end

affect! = function (integrator)
  integrator.u = integrator.u + 2
end

rootfind = true
save_positions = (true,true)
callback = ContinuousCallback(condtion,affect!,rootfind,save_positions)

sol = solve(prob,Tsit5(),callback=callback)

f = function (t,u,du)
  du[1] = - u[1] + sin(t)
end


prob = ODEProblem(f,[1.0],(0.0,10.0))

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  t - 2.95
end

affect! = function (integrator)
  integrator.u = integrator.u + 2
end

rootfind = true
save_positions = (true,true)
callback = ContinuousCallback(condtion,affect!,rootfind,save_positions)

sol = solve(prob,Tsit5(),callback=callback,abstol=1e-8,reltol=1e-6)

#=
f = @ode_def BallBounce begin
  dy =  v
  dv = -g
end g=9.81
=#

f = function (t,u,du)
  du[1] = u[2]
  du[2] = -9.81
end

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  u[1]
end

affect! = nothing
affect_neg! = function (integrator)
  integrator.u[2] = -integrator.u[2]
end

interp_points = 10
rootfind = true
save_positions = (true,true)
callback = ContinuousCallback(condtion,affect!,rootfind,save_positions,affect_neg! = affect_neg!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = ODEProblem(f,u0,tspan)


sol = solve(prob,Tsit5(),callback=callback,adaptive=false,dt=1/4)
#plot(sol,denseplot=true)

sol = solve(prob,Vern6(),callback=callback)
#plot(sol,denseplot=true)
sol = solve(prob,BS3(),callback=callback)

sol33 = solve(prob,Vern7(),callback=callback)

bounced = ODEProblem(f,sol[8],(0.0,1.0))
sol_bounced = solve(bounced,Vern6(),callback=callback,dt=sol.t[9]-sol.t[8])
#plot(sol_bounced,denseplot=true)
sol_bounced(0.04) # Complete density
@test maximum(maximum.(map((i)->sol.k[9][i]-sol_bounced.k[2][i],1:length(sol.k[9])))) == 0


sol2= solve(prob,Vern6(),callback=callback,adaptive=false,dt=1/2^4)
#plot(sol2)

sol2= solve(prob,Vern6())

sol3= solve(prob,Vern6(),saveat=[.5])

## Saving callback

condtion = function (t,u,integrator)
  true
end
affect! = function (integrator) end
save_positions = (true,false)
saving_callback = DiscreteCallback(condtion,affect!,save_positions)

sol4 = solve(prob,Tsit5(),callback=saving_callback)

@test sol2(3) ≈ sol(3)

affect! = function (integrator)
  u_modified!(integrator,false)
end
saving_callback2 = DiscreteCallback(condtion,affect!,save_positions)
sol4 = solve(prob,Tsit5(),callback=saving_callback2)

cbs = CallbackSet(saving_callback,saving_callback2)
sol4_extra = solve(prob,Tsit5(),callback=cbs)

@test length(sol4_extra) == 2length(sol4) - 1

condtion= function (t,u,integrator)
  u[1]
end

affect! = function (integrator)
  terminate!(integrator)
end

interp_points = 10
rootfind = true
save_positions = (true,true)
terminate_callback = ContinuousCallback(condtion,affect!,rootfind,save_positions)

tspan2 = (0.0,Inf)
prob2 = ODEProblem(f,u0,tspan2)

sol5 = solve(prob2,Tsit5(),callback=terminate_callback)

@test sol5[end][1] < 2e-13
@test sol5.t[end] ≈ sqrt(50*2/9.81)

affect2! = function (integrator)
  if integrator.t > 4
    terminate!(integrator)
  else
    integrator.u[2] = -integrator.u[2]
  end
end
terminate_callback2 = ContinuousCallback(condtion,affect2!,rootfind,save_positions)


sol5 = solve(prob2,Vern7(),callback=terminate_callback2)

@test sol5[end][1] < 1e-12
@test sol5.t[end] ≈ 3*sqrt(50*2/9.81)

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  t-4
end

affect! = function (integrator)
  terminate!(integrator)
end

interp_points = 10
rootfind = true
save_positions = (true,true)
terminate_callback3 = ContinuousCallback(condtion,affect!,rootfind,save_positions)

bounce_then_exit = CallbackSet(callback,terminate_callback3)

sol6 = solve(prob2,Vern7(),callback=bounce_then_exit)

@test sol6[end][1] > 0
@test sol6.t[end] ≈ 4

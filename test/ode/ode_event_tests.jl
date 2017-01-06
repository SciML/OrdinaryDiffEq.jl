using OrdinaryDiffEq, RecursiveArrayTools, DiffEqBase, Base.Test, Roots#, ParameterizedFunctions

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

affect! = function (integrator)
  integrator.u[2] = -integrator.u[2]
end

interp_points = 10
rootfind = true
save_positions = (true,true)
callback = Callback(condtion,affect!,rootfind,interp_points,save_positions)

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

condtion = true

affect! = function (integrator) end

interp_points = 0
rootfind = false
save_positions = (true,false)
saving_callback = Callback(condtion,affect!,rootfind,interp_points,save_positions)

sol4 = solve(prob,Tsit5(),callback=saving_callback)

@test sol2(3) ≈ sol(3)


condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  u[1]
end

affect! = function (integrator)
  terminate!(integrator)
end

interp_points = 10
rootfind = true
save_positions = (true,true)
terminate_callback = Callback(condtion,affect!,rootfind,interp_points,save_positions)

tspan2 = (0.0,Inf)
prob2 = ODEProblem(f,u0,tspan2)

sol5 = solve(prob2,Tsit5(),callback=terminate_callback)

@test sol5[end][1] < 2e-13
@test sol5.t[end] ≈ sqrt(50*2/9.81)

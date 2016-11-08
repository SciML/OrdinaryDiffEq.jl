using OrdinaryDiffEq, NLsolve#, ParameterizedFunctions

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

function event_f(t,u) # Event when event_f(t,u,k) == 0
  u[1]
end

function apply_event!(u,cache)
  u[2] = -u[2]
end

const dt_safety = 1
const interp_points = 10
const terminate_on_event = false
callback = @ode_callback begin
  @ode_event event_f apply_event! true interp_points terminate_on_event dt_safety
end

u0 = [50.0,0.0]
tspan = [0;15.0]
prob = ODEProblem(f,u0,tspan)


sol = solve(prob,callback=callback)
#Plots.plotly()
#plot(sol,denseplot=true)

sol = solve(prob,Vern6,callback=callback)
#plot(sol,denseplot=true)

bounced = ODEProblem(f,sol[8],[0;1.0])
sol_bounced = solve(bounced,Vern6,callback=callback,dt=sol.t[9]-sol.t[8])
#plot(sol_bounced,denseplot=true)
sol_bounced(0.04) # Complete density
bool1 = maximum(maximum.(map((i)->sol.k[9][i]-sol_bounced.k[2][i],1:length(sol.k[9])))) == 0


sol2= solve(prob,Vern6,callback=callback,adaptive=false,dt=1/2^4)
#plot(sol2)

sol2= solve(prob,Vern6)

sol3= solve(prob,Vern6,saveat=[.5])

default_callback = @ode_callback begin
  @ode_savevalues
end

sol4 = solve(prob,callback=default_callback)

bool2 = sol2(3) ≈ sol(3)

terminate_callback = @ode_callback begin
  @ode_event event_f apply_event! true interp_points true dt_safety
end

tspan2 = [0.0;Inf]
prob2 = ODEProblem(f,u0,tspan2)

sol5 = solve(prob2,callback=terminate_callback)

bool3 = sol5[end][1] < 1e-13

bool1 && bool2 && bool3

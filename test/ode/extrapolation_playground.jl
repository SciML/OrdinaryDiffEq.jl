using Revise, OrdinaryDiffEq, DiffEqBase, Plots

u0 = 1.
u0! = [1.]

f(u,p,t) = u .* cos(t)

tspan = (0. , 4. *pi)
uAnalytic(t) = exp.(sin.(t))
prob = ODEProblem(f,u0,tspan)

f!(du,u,p,t) = begin
  du[:] = cos(t) .* u[:]
end
prob! = ODEProblem(f!,u0!,tspan)

appr = solve(prob,ExtrapolationMidpointDeuflhard())
appr! = solve(prob!,ExtrapolationMidpointDeuflhard())

plt = plot(appr.t,appr.u)
plot!(plt, t->uAnalytic(t))

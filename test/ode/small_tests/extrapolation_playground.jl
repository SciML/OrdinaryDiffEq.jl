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

apprD = solve(prob,ExtrapolationMidpointDeuflhard(), dense = false)
apprD! = solve(prob!,ExtrapolationMidpointDeuflhard(), dense = false)

apprHW = solve(prob,ExtrapolationMidpointHairerWanner(), dense = false)
apprHW! = solve(prob!,ExtrapolationMidpointHairerWanner(), dense = false)

plt = plot(apprD, label = "D")
plot!(plt,apprD!, label = "D!")
plot!(plt,apprHW, label = "HW")
plot!(plt,apprHW!, label = "HW!")
plot!(plt, t->uAnalytic(t))

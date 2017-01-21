using OrdinaryDiffEq, DiffEqBase

NON_IMPLICIT_ALGS = filter((x)->isleaftype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEqAdaptiveAlgorithm)))

const α = 0.3
f = function (t,u,du)
  for i in 1:length(u)
    du[i] = α*u[i]
  end
end

condition = function (t,u,integrator)
  1-maximum(u)
end

affect! = function (integrator)
  u = integrator.u
  resize!(integrator,length(u)+1)
  maxidx = findmax(u)[2]
  Θ = rand()
  u[maxidx] = Θ
  u[end] = 1-Θ
  nothing
end

rootfind = true
save_positions = (true,true)
callback = ContinuousCallback(condition,affect!,rootfind,save_positions)

u0 = [0.2]
tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),callback=callback)

#=
using Plots
plot(sol,vars=(0,1),plotdensity=10000)

plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time")
ts = linspace(0,10,100)
plot(ts,map((x)->x[1],sol.(ts)),lw=3,
     ylabel="Amount of X in Cell 1",xlabel="Time")
=#

for alg in NON_IMPLICIT_ALGS
  if !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    println(alg)
    sol = solve(prob,alg(),callback=callback,dt=1/10)
  end
end

callback_no_interp = ContinuousCallback(condition,affect!,rootfind,save_positions)


sol = solve(prob,Tsit5(),callback=callback_no_interp,dense=false)

for alg in NON_IMPLICIT_ALGS
  if !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    println(alg)
    sol = solve(prob,alg(),callback=callback,dt=1/10)
  end
end

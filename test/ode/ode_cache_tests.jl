using OrdinaryDiffEq, DiffEqBase

NON_IMPLICIT_ALGS = filter((x)->isleaftype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEqAdaptiveAlgorithm)))

f = function (t,u,du)
  for i in 1:length(u)
    du[i] = (0.3/length(u))*u[i]
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

callback = ContinuousCallback(condition,affect!)

u0 = [0.2]
tspan = (0.0,100.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),callback=callback)

# Chunk size must be fixed since otherwise it's dependent on size
# when the size is less than 10, so errors here

sol = solve(prob,ImplicitEuler(chunk_size=1),callback=callback,dt=1/10)

sol = solve(prob,Trapezoid(chunk_size=1),callback=callback,dt=1/10)

sol = solve(prob,Rosenbrock23(chunk_size=1),callback=callback,dt=1/10)

sol = solve(prob,Rosenbrock32(chunk_size=1),callback=callback,dt=1/10)


#=
using Plots; pyplot()
p1 = plot(sol,vars=(0,1),plotdensity=10000,title="Amount of X in Cell 1")
scatter!(sol,denseplot=false)
p2 = plot(sol.t,map((x)->length(x),sol[:]),lw=3,
     ylabel="Number of Cells",xlabel="Time")
plot(p1,p2,layout=(2,1),size=(600,1000))
savefig("cell.pdf")
=#

for alg in NON_IMPLICIT_ALGS
  println(alg)
  sol = solve(prob,alg(),callback=callback,dt=1/10)
end

#=
sol = solve(prob,Tsit5(),callback=callback_no_interp,dense=false)

for alg in NON_IMPLICIT_ALGS
  println(alg)
  sol = solve(prob,alg(),callback=callback,dt=1/10)
end
=#

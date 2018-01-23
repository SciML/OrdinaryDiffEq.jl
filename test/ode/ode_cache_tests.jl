using OrdinaryDiffEq, DiffEqBase

NON_IMPLICIT_ALGS = filter((x)->isleaftype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEq.OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm)))

f = function (du,u,p,t)
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
tspan = (0.0,50.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),callback=callback)

# Chunk size must be fixed since otherwise it's dependent on size
# when the size is less than 10, so errors here

sol = solve(prob,GenericImplicitEuler(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(chunk_size=1)),callback=callback,dt=1/2)

sol = solve(prob,GenericTrapezoid(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(chunk_size=1)),callback=callback,dt=1/2)

sol = solve(prob,Rosenbrock23(chunk_size=1),callback=callback,dt=1/2)

sol = solve(prob,Rosenbrock32(chunk_size=1),callback=callback,dt=1/2)

for alg in CACHE_TEST_ALGS
  sol = solve(prob,alg,callback=callback,dt=1/2)
end

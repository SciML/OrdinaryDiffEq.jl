using OrdinaryDiffEq, DiffEqBase, DiffEqCallbacks, Test
using Random
Random.seed!(213)
CACHE_TEST_ALGS = [Euler(),Midpoint(),RK4(),SSPRK22(),SSPRK33(),
  ORK256(), CarpenterKennedy2N54(), HSLDDRK64(), DGLDDRK73_C(),
  CFRLDDRK64(), TSLDDRK74(),
  CKLLSRK43_2(),
  ParsaniKetchesonDeconinck3S32(),
  BS3(),BS5(),DP5(),DP8(),Feagin10(),Feagin12(),Feagin14(),TanYam7(),
  Tsit5(),TsitPap8(),Vern6(),Vern7(),Vern8(),Vern9(),OwrenZen3(),OwrenZen4(),OwrenZen5()]

using InteractiveUtils

NON_IMPLICIT_ALGS = filter((x)->isconcretetype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEq.OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm)))

f = function (du,u,p,t)
  for i in 1:length(u)
    du[i] = (0.3/length(u))*u[i]
  end
end

condition = function (u,t,integrator)
  1-maximum(u)
end

affect! = function (integrator)
  u = integrator.u
  resize!(integrator,length(u)+1)
  maxidx = findmax(u)[2]
  Θ = rand()/5 + 0.25
  u[maxidx] = Θ
  u[end] = 1-Θ
  nothing
end

callback = ContinuousCallback(condition,affect!)

u0 = [0.2]
tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan)

println("Check for stochastic errors")
for i in 1:10
  @test_nowarn sol = solve(prob,Tsit5(),callback=callback)
end

println("Check some other integrators")
@test_nowarn sol = solve(prob,GenericImplicitEuler(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(chunk_size=1)),callback=callback,dt=1/2)
@test_nowarn sol = solve(prob,GenericTrapezoid(nlsolve=OrdinaryDiffEq.NLSOLVEJL_SETUP(chunk_size=1)),callback=callback,dt=1/2)
@test_nowarn sol = solve(prob,Rosenbrock23(chunk_size=1),callback=callback,dt=1/2)
@test_nowarn sol = solve(prob,Rosenbrock32(chunk_size=1),callback=callback,dt=1/2)

for alg in CACHE_TEST_ALGS
  @test_nowarn sol = solve(prob,alg,callback=callback,dt=1/2)
end

@test_nowarn sol = solve(prob,Rodas4(chunk_size=1),callback=callback,dt=1/2)
@test_nowarn sol = solve(prob,Rodas5(chunk_size=1),callback=callback,dt=1/2)

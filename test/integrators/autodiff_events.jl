using DiffEqSensitivity
using OrdinaryDiffEq, Calculus, Test

function f(du,u,p,t)
  du[1] = -p[1]
  du[2] = p[2]
end

cb = ContinuousCallback((u,t,i) -> u[1], (integrator)->(println("Stopped.");integrator.p[2]=0.0))
p = [2.0, 1.0]

function test_f(p)
  prob = ODEProblem(f,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
  solve(prob,Tsit5(),abstol=1e-14,reltol=1e-14,callback=cb,save_everystep=false)[end]
end
p = [2.0, 1.0]
findiff = Calculus.finite_difference_jacobian(test_f,p)
findiff

using ForwardDiff
ad = ForwardDiff.jacobian(test_f,p)
ad

@test ad ≈ findiff

using ParameterizedFunctions, Test
using OrdinaryDiffEq, Calculus, ForwardDiff

f = @ode_def begin
    dx = -a
    dy = b
end a b

cb = ContinuousCallback((u,t,i) -> u[1], (integrator)->(println("Stopped.");integrator.p[2]=zero(integrator.p[2])))
function test_f(p)
  prob = ODEProblem(f,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
  integrator = init(prob,Tsit5(),abstol=1e-14,reltol=1e-14,callback=cb)
  step!(integrator)
  solve!(integrator).u[end]
end
p = [2.0, 1.0]
findiff = Calculus.finite_difference_jacobian(test_f,p)
fordiff = ForwardDiff.jacobian(test_f,p)
@test findiff ≈ fordiff

f2 = @ode_def begin
    dx = -x
    dy = b
end a b

function test_f2(p)
  prob = ODEProblem(f2,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
  integrator = init(prob,Tsit5(),abstol=1e-14,reltol=1e-14,callback=cb)
  step!(integrator)
  solve!(integrator).u[end]
end
p = [2.0, 1.0]
findiff = Calculus.finite_difference_jacobian(test_f2,p)
fordiff = ForwardDiff.jacobian(test_f2,p)
@test findiff ≈ fordiff

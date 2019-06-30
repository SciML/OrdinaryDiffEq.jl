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

# Gradients and Hessians

function myobj(θ)
  f(u,p,t) = -θ[1]*u
  u0, _ = promote(10.0, θ[1])
  prob = ODEProblem(f, u0, (0.0, 1.0))
  sol = solve(prob, Tsit5())
  diff = sol.u - 10*exp.(-sol.t)
  return diff'diff
end

ForwardDiff.gradient(myobj, [1.0])
ForwardDiff.hessian(myobj, [1.0])

function myobj2(θ)
  f(du,u,p,t) = (du[1]=-θ[1]*u[1])
  u0, _ = promote(10.0, θ[1])
  prob = ODEProblem(f, [u0], (0.0, 1.0))
  sol = solve(prob, Tsit5())
  diff = sol[:,1] .- 10 .*exp.(-sol.t)
  return diff'diff
end

ForwardDiff.gradient(myobj2, [1.0])
ForwardDiff.hessian(myobj2, [1.0])

function myobj3(θ)
  f(u,p,t) = -θ[1]*u
  u0, _ = promote(10.0, θ[1])
  tspan_end, _ = promote(1.0, θ[1])
  prob = ODEProblem(f, u0, (0.0, tspan_end))
  sol = solve(prob, Tsit5())
  diff = sol.u - 10*exp.(-sol.t)
  return diff'diff
end

ForwardDiff.gradient(myobj3, [1.0])
ForwardDiff.hessian(myobj3, [1.0])

function myobj4(θ)
  f(du,u,p,t) = (du[1] = -θ[1]*u[1])
  u0, _ = promote(10.0, θ[1])
  tspan_end, _ = promote(1.0, θ[1])
  prob = ODEProblem(f, [u0], (0.0, tspan_end))
  sol = solve(prob, Tsit5())
  diff = sol[:,1] .- 10 .* exp.(-sol.t)
  return diff'diff
end

ForwardDiff.gradient(myobj4, [1.0])
ForwardDiff.hessian(myobj4, [1.0])

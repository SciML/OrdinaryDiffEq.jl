using Test
using OrdinaryDiffEq, Calculus, ForwardDiff

function f(du,u,p,t)
  du[1] = -p[1]
  du[2] = p[2]
end

for x in 0:0.001:5
  called = false
  function test_f(p)
    cb = ContinuousCallback((u,t,i) -> u[1], (integrator)->(called=true;integrator.p[2]=zero(integrator.p[2])))
    prob = ODEProblem(f,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
    integrator = init(prob,Tsit5(),abstol=1e-14,reltol=1e-14,callback=cb)
    step!(integrator)
    solve!(integrator).u[end]
  end
  p = [2.0, x]
  called = false
  findiff = Calculus.finite_difference_jacobian(test_f,p)
  @test called
  called = false
  fordiff = ForwardDiff.jacobian(test_f,p)
  @test called
  @test findiff ≈ fordiff
end

function f2(du,u,p,t)
  du[1] = -u[2]
  du[2] = p[2]
end

for x in 2.1:0.001:5
  called = false
  function test_f2(p)
    cb = ContinuousCallback((u,t,i) -> u[1], (integrator)->(called=true;integrator.p[2]=zero(integrator.p[2])))
    prob = ODEProblem(f2,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
    integrator = init(prob,Tsit5(),abstol=1e-12,reltol=1e-12,callback=cb)
    step!(integrator)
    solve!(integrator).u[end]
  end
  p = [2.0, x]
  findiff = Calculus.finite_difference_jacobian(test_f2,p)
  @test called
  called = false
  fordiff = ForwardDiff.jacobian(test_f2,p)
  @test called
  @test findiff ≈ fordiff
end

#=
#x = 2.0 is an interesting case

x = 2.0

function test_f2(p)
  cb = ContinuousCallback((u,t,i) -> u[1], (integrator)->(@show(x,integrator.t);called=true;integrator.p[2]=zero(integrator.p[2])))
  prob = ODEProblem(f2,eltype(p).([1.0,0.0]),eltype(p).((0.0,1.0)),copy(p))
  integrator = init(prob,Tsit5(),abstol=1e-12,reltol=1e-12,callback=cb)
  step!(integrator)
  solve!(integrator).u[end]
end

p = [2.0, x]
findiff = Calculus.finite_difference_jacobian(test_f2,p)
@test called
called = false
fordiff = ForwardDiff.jacobian(test_f2,p)
@test called

# At that value, it shouldn't be called, but a small perturbation will make it called, so finite difference is wrong!
=#

for x in 1.0:0.001:2.5
  function lotka_volterra(du,u,p,t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α*x - β*x*y
    du[2] = dy = -δ*y + γ*x*y
  end
  u0 = [1.0,1.0]
  tspan = (0.0,10.0)
  p = [x,1.0,3.0,1.0]
  prob = ODEProblem(lotka_volterra,u0,tspan,p)
  sol = solve(prob,Tsit5())

  called=false
  function test_lotka(p)
    cb = ContinuousCallback((u,t,i) -> u[1]-2.5, (integrator)->(called=true;integrator.p[4]=1.5))
    prob = ODEProblem(lotka_volterra,eltype(p).([1.0,1.0]),eltype(p).((0.0,10.0)),copy(p))
    integrator = init(prob,Tsit5(),abstol=1e-12,reltol=1e-12,callback=cb)
    step!(integrator)
    solve!(integrator).u[end]
  end

  findiff = Calculus.finite_difference_jacobian(test_lotka,p)
  @test called
  called = false
  fordiff = ForwardDiff.jacobian(test_lotka,p)
  @test called
  @test findiff ≈ fordiff
end

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

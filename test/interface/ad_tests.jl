using Test
using OrdinaryDiffEq, Calculus, ForwardDiff

function f(du,u,p,t)
  du[1] = -p[1]
  du[2] = p[2]
end

for x in 0:0.001:5
  called = false
  if x in [1.0, 2.0, 3.0, 4.0, 5.0]
    print("AD Ping $x")
  end
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
  if x in [3.0, 4.0, 5.0]
    print("AD Ping $x")
  end
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
  if x in [1.5, 2.0, 2.5]
    print("AD Ping $x")
  end
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


f1s = function (du,u,p,t)
    du[1] = p[2] * u[1]
    du[2] = p[3] * u[2]
    du[3] = p[4] * u[3]
    du[4] = p[5] * u[4]
    nothing
end
f2s = function (du,u,p,t)
    du[1] = p[1]*u[2]
    du[2] = p[1]*u[3]
    du[3] = p[1]*u[4]
    du[4] = p[1]*u[1]
    nothing
end
u0 = [3.4, 3.3, 3.2, 3.1]
params = [0.002, -0.005, -0.004, -0.003, -0.002]
tspan = (7.0, 84.0)
times = collect(minimum(tspan):0.5:maximum(tspan))
prob = SplitODEProblem(f1s,f2s, u0, tspan, params)
sol2 = solve(prob, KenCarp4(); dt=0.5, saveat=times)

function difffunc(p)
    tmp_prob = remake(prob,p=p)
    vec(solve(tmp_prob,KenCarp4(),saveat=times))
end
ForwardDiff.jacobian(difffunc,ones(5))

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1221

f_a  = function (du, u, p, t)
    du[1] = -p[1]*u[1] + exp(-t)
end

of_a = p -> begin
    u0  = [0.0]
    tspan = (0.0, 5.0)
    prob  = ODEProblem(f_a, u0, tspan, p)
    # sol = solve(prob, Tsit5())                      # works
    # sol = solve(prob, Rodas5(autodiff=false))       # works
    sol = solve(prob, Rodas5(autodiff=true),abstol=1e-14,reltol=1e-14)          # fails
    return sum(t -> abs2(t[1]), sol([1.0, 2.0, 3.0]))
end

@test !iszero(ForwardDiff.gradient(t -> of_a(t), [1.0]))
@test ForwardDiff.gradient(t -> of_a(t), [1.0]) ≈ FiniteDiff.finite_difference_gradient(t -> of_a(t), [1.0]) rtol=1e-5

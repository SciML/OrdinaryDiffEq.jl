using Test
using OrdinaryDiffEq, Calculus, ForwardDiff, FiniteDiff, LinearAlgebra

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

# From https://discourse.julialang.org/t/issue-with-callbacks-when-solving-differential-equations-having-dual-number-states/41883/6

#                     val    d/dy0  d/dvy0 d/dgy
y  = ForwardDiff.Dual(0.6,   1.0,   0.0,   0.0  )
vy = ForwardDiff.Dual(0.0,   0.0,   1.0,   0.0  )
gy = ForwardDiff.Dual(0.005, 0.0,   0.0,   1.0  )
u0 = [y;vy]
tspan = (0.0, 2.0)
function st!(du, u, p, t)
    @inbounds begin
        du[1] = u[2]
        du[2] = (1-u[1]) - p
    end
    return nothing
end
cnd(u, t, integ) = u[1] - 1
CC = ContinuousCallback(cnd, terminate!, nothing, abstol=0, reltol=0)

prb = ODEProblem(st!, u0, tspan, gy)
sl1 = solve(prb, Vern9(), abstol=1e-12, reltol=1e-12)
sl2 = solve(prb, Vern9(), abstol=1e-12, reltol=1e-12, callback=CC)

s1 = sl1(sl2.t[end].value)
s2 = sl2[end]
@test s1 ≈ s2

s3 = sl1(sl2.t[end-1].value)
s4 = sl2[end-1]
@test s3 ≈ s4

s5 = sl1(sl2.t[end-2].value)
s6 = sl2[end-2]
@test s5 ≈ s6

s7 = sl1(sl2.t[end-3].value)
s8 = sl2[end-3]
@test s7 ≈ s8

s9 = sl1(sl2.t[end-4].value)
s10= sl2[end-4]
@test s9 ≈ s10

x = [0.6,0.0,0.005]
function get_endminusidx_cb(x;idx=0)
    y  = x[1]
    vy = x[2]
    gy = x[3]
    u0 = [y;vy]
    tspan = (0.0, 2.0)
    function st!(du, u, p, t)
        @inbounds begin
            du[1] = u[2]
            du[2] = (1-u[1]) - p
        end
        return nothing
    end
    cnd(u, t, integ) = u[1] - 1
    CC = ContinuousCallback(cnd, terminate!, nothing, abstol=0, reltol=0)

    prb = ODEProblem(st!, u0, tspan, gy)
    sol = solve(prb, Vern9(), abstol=1e-12, reltol=1e-12, callback = CC, saveat= ForwardDiff.value.(sl2.t))
    sol[end-idx]
end
function get_endminusidx(x;idx=0)
    y  = x[1]
    vy = x[2]
    gy = x[3]
    u0 = [y;vy]
    tspan = (0.0, 2.0)
    function st!(du, u, p, t)
        @inbounds begin
            du[1] = u[2]
            du[2] = (1-u[1]) - p
        end
        return nothing
    end
    cnd(u, t, integ) = u[1] - 1
    CC = ContinuousCallback(cnd, terminate!, nothing, abstol=0, reltol=0)

    prb = ODEProblem(st!, u0, tspan, gy)
    sol = solve(prb, Vern9(), abstol=1e-12, reltol=1e-12, saveat= ForwardDiff.value.(sl2.t))
    sol[end-idx]
end
# The first two need to use the callback, since that makes finite difference
# know that u[1] does not change at the end
enddiffs = FiniteDiff.finite_difference_jacobian(get_endminusidx_cb,x)
@test enddiffs ≈ reduce(hcat,Array.(ForwardDiff.partials.(s2)))' atol=1e-7

enddiffs = FiniteDiff.finite_difference_jacobian(x->get_endminusidx_cb(x,idx=1),x)
@test enddiffs ≈ reduce(hcat,Array.(ForwardDiff.partials.(s4)))' atol=1e-7

# These now use the no-callback version since that seems to be more stable
# for finite differencing
enddiffs = FiniteDiff.finite_difference_jacobian(x->get_endminusidx(x,idx=2),x)
@test enddiffs ≈ reduce(hcat,Array.(ForwardDiff.partials.(s6)))' atol=1e-5

enddiffs = FiniteDiff.finite_difference_jacobian(x->get_endminusidx(x,idx=3),x)
@test enddiffs ≈ reduce(hcat,Array.(ForwardDiff.partials.(s8)))' atol=1e-5

enddiffs = FiniteDiff.finite_difference_jacobian(x->get_endminusidx(x,idx=4),x)
@test enddiffs ≈ reduce(hcat,Array.(ForwardDiff.partials.(s10)))' atol=1e-5

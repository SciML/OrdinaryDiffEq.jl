using DiffEqBase, OrdinaryDiffEq, Base.Test

type SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end

f = function (t,u,du)
    du[1] = -0.5*u[1] + u.f1
    du[2] = -0.5*u[2]
end

const tstop1 = [5.]
const tstop2 = [8.]
const tstop = [5.;8.]

condition = function (t,u,integrator)
  t in tstop1
end

affect! = function (integrator)
  for c in user_cache(integrator)
    c.f1 = 1.5
  end
end

save_positions = (true,true)

cb = DiscreteCallback(condition, affect!; save_positions=save_positions)

condition2 = function (t,u,integrator)
  t in tstop2
end

affect2! = function (integrator)
  for c in user_cache(integrator)
    c.f1 = -1.5
  end
end

save_positions = (true,true)

cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)

cbs = CallbackSet(cb,cb2)

u0 = SimType{Float64}([10;10], 0.0)
prob = ODEProblem(f,u0,(0.0,10.0))
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)

sol(1.5:0.5:2.5)

@test [sol[i].f1 for i in eachindex(sol)] == [zeros(9);1.5*ones(5);-1.5*ones(4)]

using DiffEqBase, OrdinaryDiffEq

A = diagm([0.3,0.6,0.9])
B = [1 2 3].'
C = [1/3 1/3 1/3]

type SimType2{T} <: DEDataVector{T}
    x::Vector{T}
    y::Vector{T}
    u::Vector{T}
end

function mysystem(t,x,dx,u)
    ucalc = u(t,x)
    x.u = ucalc
    x.y = C*x.x
    dx[:] = A*x.x + B*x.u
end

input = (t,x)->(1*one(t)≤t≤2*one(t)?[one(t)]:[zero(t)])
prob = DiscreteProblem((t,x,dx)->mysystem(t,x,dx,input), SimType2(zeros(3), zeros(1), zeros(1)), (0//1,4//1))
sln = solve(prob, FunctionMap(scale_by_time=false), dt = 1//10)

u1 = [sln[idx].u for idx in 1:length(sln)]
u2 = [sln(t).u for t in linspace(0,4,41)]
@test any(x->x[1]>0, u1)
@test any(x->x[1]>0, u2)

sln = solve(prob, FunctionMap(scale_by_time=true), dt = 1//10)

u1 = [sln[idx].u for idx in 1:length(sln)]
u2 = [sln(t).u for t in linspace(0,4,41)]
@test any(x->x[1]>0, u1)
@test any(x->x[1]>0, u2)

sln = solve(prob, Euler(), dt = 1//10)

@test u1 == [sln[idx].u for idx in 1:length(sln)] # Show that discrete is the same
u1 = [sln[idx].u for idx in 1:length(sln)]
u2 = [sln(t).u for t in linspace(0,4,41)]
@test any(x->x[1]>0, u1)
@test any(x->x[1]>0, u2)

sln = solve(prob, DP5(), dt = 1//10, adaptive=false)

u1 = [sln[idx].u for idx in 1:length(sln)]
u2 = [sln(t).u for t in linspace(0,4,41)]
@test any(x->x[1]>0, u1)
@test any(x->x[1]>0, u2)

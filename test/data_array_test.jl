using OrdinaryDiffEq, Base.Test

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

using OrdinaryDiffEq

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





######################
# DEDataMatrix

type SimTypeg{T,T2} <: DEDataMatrix{T}
    x::Array{T,2} # two dimensional
    f1::T2
end

const tstop1 = [10.0]
const tstop2 = [300.]

function mat_condition(t,u,integrator)
  t in tstop1
end

function mat_condition2(t,u,integrator)
  t in tstop2
end

function mat_affect!(integrator)
  for c in user_cache(integrator)
    c.f1 = +1.0
  end
#  integrator.u[1,1] = 0.001
end

function mat_affect2!(integrator)
  for c in user_cache(integrator)
    c.f1 = 0.0
  end
end

save_positions = (true,true)
cb = DiscreteCallback(mat_condition, mat_affect!, save_positions=save_positions)
save_positions = (false,true)
cb2 = DiscreteCallback(mat_condition2, mat_affect2!, save_positions=save_positions)
cbs = CallbackSet(cb,cb2)


function sigmoid(t,u,du)
  du[1,1] = 0.01*u[1,1]*(1-u[1,1]/20)
  du[1,2] = 0.01*u[1,2]*(1-u[1,2]/20)
  du[2,1] = 0.01*u[2,1]*(1-u[2,1]/20)
  du[2,2] = u.f1*du[1,1]
end

u0 = SimTypeg(fill(0.00001,2,2),0.0)
tspan = (0.0,3000.0)
prob = ODEProblem(sigmoid,u0,tspan)

const tstop =[tstop1;tstop2]
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)
sol = solve(prob,Rodas4(),callback = cbs, tstops=tstop)
sol = solve(prob,Kvaerno3(),callback = cbs, tstops=tstop)
@test_broken sol = solve(prob,Rodas4(autodiff=false),callback = cbs, tstops=tstop)
@test_broken sol = solve(prob,Kvaerno3(autodiff=false),callback = cbs, tstops=tstop)

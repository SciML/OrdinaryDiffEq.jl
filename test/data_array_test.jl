using DiffEqBase, OrdinaryDiffEq

type SimType{T} <: DEDataArray{T}
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

cb = DiscreteCallback(condition, affect!, save_positions)

condition2 = function (t,u,integrator)
  t in tstop2
end

affect2! = function (integrator)
  for c in user_cache(integrator)
    c.f1 = -1.5
  end
end

save_positions = (false,true)

cb2 = DiscreteCallback(condition2, affect2!, save_positions)

cbs = CallbackSet(cb,cb2)

u0 = SimType{Float64}([10;10], 0.0)
prob = ODEProblem(f,u0,(0.0,10.0))
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)

sol(1.5:0.5:2.5)

@test [sol[i].f1 for i in eachindex(sol)] == [zeros(9);1.5*ones(5);-1.5*ones(4)]

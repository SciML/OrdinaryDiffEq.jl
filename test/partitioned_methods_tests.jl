using DiffEqBase, OrdinaryDiffEq, Base.Test, RecursiveArrayTools

u0 = 3ones(4)
v0 = ones(4)
f1 = function (t,u,v,du)
  du .= v
end
f2 = function (t,u,v,dv)
  dv .= -2u
end

prob = PartitionedODEProblem((f1,f2),(u0,v0),(0.0,5.0))

sol = solve(prob,SymplecticEuler(),dt=1/100)

interp_time = 0:0.001:5
interps = sol(interp_time)

prob = SecondOrderODEProblem(f2,u0,v0,(0.0,5.0))
sol2 = solve(prob,SymplecticEuler(),dt=1/100)

@test sol[end][1] == sol2[end][1]
@test sol[end][5] == sol2[end][5]

f = function (t,u,du)
  du.x[1] .= u.x[2]
  du.x[2] .= -2u.x[1]
end

u = ArrayPartition((u0,v0))

prob = ODEProblem(f,u,(0.0,5.0))

sol = solve(prob,Euler(),dt=1/100)

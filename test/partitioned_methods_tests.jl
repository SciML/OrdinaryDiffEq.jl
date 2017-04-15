using DiffEqBase, OrdinaryDiffEq

u0 = 3ones(4)
v0 = ones(4)
f1 = function (t,u,v,du)
  du .= v
end
f2 = function (t,u,v,du)
  dv .= -2u
end

prob = PartitionedODEProblem((f1,f2),(u0,v0),(0.0,1.0))

integrator = init(prob,Discrete())

sol = solve(prob,Discrete())

using OrdinaryDiffEq, DiffEqProblemLibrary
prob = prob_ode_linear

sol = solve(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
sol(0.9)

integrator = init(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)

for i in integrator
  @show i.t,i.u
end
integrator.sol(0.9)

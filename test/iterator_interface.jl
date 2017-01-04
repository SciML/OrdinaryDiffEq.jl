using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test
prob = prob_ode_linear

sol = solve(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
sol(0.9)

integrator = init(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
start(integrator)
step(integrator)
solve!(integrator)

@test integrator.sol(0.9) == sol(0.9)

integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5],advance_to_tstop=true)
for i in integrator
  @test i.t ∈ [0.5,1.0]
end

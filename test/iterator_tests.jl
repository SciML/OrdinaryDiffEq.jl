using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test
prob = prob_ode_linear

sol = solve(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
sol(0.9)

integrator = init(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
start(integrator)
step(integrator)
@test integrator.iter == 1
solve!(integrator)
@test integrator.t == 1.0

push!(integrator.opts.tstops,5.0)
integrator.opts.advance_to_tstop=true
step(integrator)
@test integrator.t == 5.0
integrator.opts.advance_to_tstop=false
step(integrator)
@test integrator.t > 5

@test integrator.sol(0.9) == sol(0.9)

integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5],advance_to_tstop=true)
for i in integrator
  @test i.t ∈ [0.5,1.0]
end

integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5],advance_to_tstop=true,stop_at_next_tstop=true)
for i in integrator
  @test i.t == 0.5
end

#=
using Plots
integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
pyplot(show=true)
plot(integrator)
for i in integrator
  plot!(integrator,vars=(0,1),legend=false)
end
step(integrator); plot!(integrator,vars=(0,1),legend=false)
=#

using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test, DiffEqBase
prob = prob_ode_linear

sol = solve(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
sol(0.9)

integrator = init(prob,BS3();dt=1//2^(4),tstops=[0.5],saveat=0:0.01:1)
step!(integrator)
@test integrator.iter == 1
solve!(integrator)
@test integrator.t == 1.0
integrator(0.95)
integrator.tprev
integrator.t

push!(integrator.opts.tstops,5.0)
integrator.opts.advance_to_tstop=true
step!(integrator)
@test integrator.t == 5.0
integrator.opts.advance_to_tstop=false
step!(integrator)
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
integrator([1.0;2.0])

integrator = init(prob,RK4();dt=1//2^(9))
for i in take(integrator,12) end
@test integrator.iter == 12
for i in take(integrator,12) end
@test integrator.iter == 24

#=
using Plots
integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
pyplot(show=true)
plot(integrator)
for i in integrator
  display(plot!(integrator,vars=(0,1),legend=false))
end
step!(integrator); plot!(integrator,vars=(0,1),legend=false)
savefig("iteratorplot.png")

integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
plot(integrator)
for i in integrator
  display(plot!(integrator,vars=(0,1),legend=false,denseplot=false))
end
step!(integrator); plot!(integrator,vars=(0,1),legend=false,denseplot=false)
=#

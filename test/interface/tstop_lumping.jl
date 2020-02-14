using OrdinaryDiffEq, Test

function affect!(integrator)
    integrator.u += 0.05
    nothing
end
cb = PeriodicCallback(affect!, 0.1; initial_affect = false)
prob = ODEProblem((u,u0,t)->1.01u,0.5,(0.0,1.0))
sim = solve(prob,Tsit5(),callback = cb)
@test sim.t[end] == 1.0
@test sim.t[end] - sim.t[end-1] > 10eps()

prob = ODEProblem((u,u0,t)->1.01u,0.5,(0.0,1.0))
sim = solve(prob,Tsit5(),callback = cb,tstops = 0.0:0.1:1.0)
@test sim.t[end] == 1.0
@test sim.t[end] - sim.t[end-1] > 10eps()

prob = ODEProblem((u,u0,t)->1.01u,0.5,(0.0,2.0))
sim = solve(prob,Tsit5(),callback = cb)
@test sim.t[end] == 2.0
@test sim.t[end] - sim.t[end-1] > 10eps()

function affect!(integrator)
    nothing
end
cb = PeriodicCallback(affect!, 0.1; initial_affect = false)
prob = ODEProblem((u,u0,t)->1.01u,0.5,(0.0,1.0))
sim = solve(prob,Tsit5())
@test sim.t[end] == 1.0
@test sim.t[end] - sim.t[end-1] > 10eps()

sim = solve(prob,Tsit5(),callback = cb)
@test sim.t[end] == 1.0
@test sim.t[end] - sim.t[end-1] > 10eps()

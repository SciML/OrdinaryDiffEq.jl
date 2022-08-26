using OrdinaryDiffEq, DiffEqCallbacks, Test
import ODEProblemLibrary: prob_ode_2Dlinear

prob = prob_ode_2Dlinear

integrator = init(prob, Tsit5(); dt = 1 // 2^(4))
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

# Implicitly test if dt resets
integrator = init(prob, Tsit5())
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

integrator = init(prob, Tsit5(); dt = 1 // 2^(4), tstops = [1 // 2], saveat = [1 // 4])
dt = integrator.dt
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

#callback test
g(u, p, t) = 2.0 * t - 2.0
u0 = 0.0
tspan = (0.0, 4.0)
prob = ODEProblem(g, u0, tspan)
saved_values = SavedValues(Float64, Float64)
cb = SavingCallback((u, t, integrator) -> u, saved_values)
integrator = init(prob, Tsit5(); dt = 1 // 2^(4), callback = cb)
dt = integrator.dt
solve!(integrator)

u = saved_values.saveval
t = saved_values.t
resize!(saved_values.t, 0)
resize!(saved_values.saveval, 0)
reinit!(integrator)
integrator.dt = dt
solve!(integrator)

@test u == saved_values.saveval
@test t == saved_values.t

@testset "set u0" begin
    prob = prob_ode_2Dlinear
    integrator = init(prob, Tsit5())
    u0 = prob.u0 .+ 1  # just make it different
    @test u0 != prob.u0
    reinit!(integrator, u0)
    @test integrator.u == u0
    @test integrator.sol.u[1] == u0
    @test integrator.sol.interp.timeseries[1] == u0
end

@testset "set u0 with save_idxs" begin
    save_idxs = [1]
    prob = prob_ode_2Dlinear
    integrator = init(prob, Tsit5(); save_idxs = save_idxs)
    u0 = prob.u0 .+ 1  # just make it different
    @test u0 != prob.u0
    reinit!(integrator, u0)
    @test integrator.u == u0
    @test integrator.sol.u[1] == u0[save_idxs]
    @test integrator.sol.interp.timeseries[1] == u0[save_idxs]
end

@testset "set t0" begin
    prob = prob_ode_2Dlinear
    integrator = init(prob, Tsit5())
    t0 = prob.tspan[1] - 1  # just make it different
    @test t0 != prob.tspan[1]
    reinit!(integrator; t0 = t0)
    @test integrator.t == t0
    @test integrator.sol.t[1] == t0
    @test integrator.sol.interp.ts[1] == t0
end

using DelayDiffEq, DDEProblemLibrary, DiffEqDevTools, DiffEqCallbacks
using OrdinaryDiffEqTsit5
using Statistics
using Test

const prob = prob_dde_constant_1delay_scalar
const alg = MethodOfSteps(Tsit5(); constrained = false)

# continuous callback
@testset "continuous" begin
    cb = ContinuousCallback(
        (u, t, integrator) -> t - 2.6, # Event when event_f(t,u,k) == 0
        integrator -> (integrator.u = -integrator.u)
    )

    sol1 = solve(prob, alg, callback = cb)
    ts = findall(x -> x ≈ 2.6, sol1.t)
    @test length(ts) == 2
    @test sol1.u[ts[1]] == -sol1.u[ts[2]]
    @test sol1(
        2.6;
        continuity = :right
    ) ≈ -sol1(2.6; continuity = :left) atol = 2.0e-5

    # fails on 32bit?!
    # see https://github.com/SciML/DelayDiffEq.jl/pull/180
    if Int === Int64
        sol2 = solve(prob, alg, callback = cb, dtmax = 0.01)
        ts = findall(x -> x ≈ 2.6, sol2.t)
        @test length(ts) == 2
        @test sol2.u[ts[1]] == -sol2.u[ts[2]]
        @test sol2(2.6; continuity = :right) ≈
            -sol2(2.6; continuity = :left)

        sol3 = appxtrue(sol1, sol2)
        @test sol3.errors[:L2] < 1.5e-2
        @test sol3.errors[:L∞] < 3.5e-2
    end
end

# discrete callback
@testset "discrete" begin
    # Automatic absolute tolerances
    cb = AutoAbstol()

    sol1 = solve(prob, alg, callback = cb)
    sol2 = solve(prob, alg, callback = cb, dtmax = 0.01)
    sol3 = appxtrue(sol1, sol2)

    @test sol3.errors[:L2] < 1.4e-3
    @test sol3.errors[:L∞] < 4.1e-3

    # Terminate early
    cb = DiscreteCallback((u, t, integrator) -> t == 4, terminate!)
    sol = @test_logs solve(prob, alg; callback = cb, tstops = [4.0])
    @test sol.t[end] == 4
end

@testset "preset callback at constant-lag propagated tstop" begin
    function lif2_Ndelay!(du, u, h, p, t)
        u1, u2 = u
        gL1, EL1, C1, _, I1, gL2, EL2, C2, _, I2, λ, T... = p

        H = Float64[h(p, t - τ)[2] for τ in T]
        append!(H, u2)

        du[1] = (-gL1 * (u1 - EL1) + I1 - λ * (u1 - mean(H))) / C1
        du[2] = (-gL2 * (u2 - EL2) + I2) / C2
    end

    fired = Float64[]
    up_Iex2!(integrator) = (push!(fired, integrator.t); integrator.p[10] = 210.0)
    down_Iex2!(integrator) = integrator.p[10] = 0.0

    cb = CallbackSet(
        ContinuousCallback(
            (u, t, integrator) -> u[1] - integrator.p[4],
            integrator -> (integrator.u[1] = integrator.p[2])
        ),
        ContinuousCallback(
            (u, t, integrator) -> u[2] - integrator.p[9],
            integrator -> (integrator.u[2] = integrator.p[7])
        ),
        PresetTimeCallback([2.0, 9.0], up_Iex2!),
        PresetTimeCallback([3.0, 10.0], down_Iex2!)
    )

    u0 = [-75.0, -75.0]
    tspan = (0.0, 12.0)
    p1 = [10.0, -75.0, 5.0, -55.0, 100.0]
    p2 = [10.0, -75.0, 5.0, -55.0, 0.0]
    delays = [0.2, 4.0]
    p = [p1; p2; [2.0]; delays]
    h(p, t) = u0

    prob = DDEProblem(lif2_Ndelay!, u0, h, tspan, p; callback = cb, constant_lags = delays)
    sol = solve(prob, MethodOfSteps(Tsit5()))

    @test sol.retcode == ReturnCode.Success
    @test any(==(9.0), fired)
end

@testset "save discontinuity" begin
    f(du, u, h, p, t) = (du .= 0)
    prob = DDEProblem(f, [0.0], nothing, (0.0, 1.0))

    condition(u, t, integrator) = t == 0.5
    global iter
    function affect!(integrator)
        integrator.u[1] += 100
        global iter = integrator.iter
    end
    cb = DiscreteCallback(condition, affect!)
    sol = solve(prob, MethodOfSteps(Tsit5()), callback = cb, tstops = [0.5])
    @test sol.t[iter + 1] == sol.t[iter + 2]
    @test sol.u[iter + 1] == [0.0]
    @test sol.u[iter + 2] != [0.0]
end

# Issue #341: PeriodicCallback support for DDEIntegrator
@testset "periodic callback" begin
    # Simple DDE with constant delay
    function f(du, u, h, p, t)
        du[1] = -u[1] + h(p, t - 1.0)[1]
    end
    h(p, t) = [1.0]
    prob = DDEProblem(f, [1.0], h, (0.0, 10.0); constant_lags = [1.0])

    # Count how many times the callback is triggered
    counter = Ref(0)
    function affect!(integrator)
        counter[] += 1
    end
    cb = PeriodicCallback(affect!, 1.0)

    sol = solve(prob, MethodOfSteps(Tsit5()), callback = cb)
    @test sol.retcode == ReturnCode.Success
    # Should be triggered approximately 10 times (at t=1,2,3,...,10)
    @test counter[] >= 9
end

using OrdinaryDiffEq, Test, Random, StaticArrays, DiffEqCallbacks
import ODEProblemLibrary: prob_ode_linear
Random.seed!(100)

@testset "Tstops Tests on the Interval [0, 1]" begin
    prob = prob_ode_linear

    sol = solve(prob, Tsit5(), dt = 1 // 2^(6), tstops = [1 / 2])
    @test 1 // 2 ∈ sol.t

    sol = solve(prob, RK4(), dt = 1 // 3, tstops = [1 / 2], adaptive = false)
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    sol = solve(
        prob, RK4(), dt = 1 // 3, tstops = [1 / 2],
        d_discontinuities = [-1 / 2, 1 / 2, 3 / 2], adaptive = false
    )
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    # TODO
    integrator = init(
        prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4],
        adaptive = false
    )

    sol = solve(prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4], adaptive = false)
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(
        prob, RK4(), tstops = [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1],
        adaptive = false
    )
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(prob, RK4(), tstops = 0:(1 // 16):1, adaptive = false)
    @test sol.t == collect(0:(1 // 16):1)

    sol = solve(prob, RK4(), tstops = range(0, stop = 1, length = 100), adaptive = false)
    @test sol.t == collect(range(0, stop = 1, length = 100))
end

@testset "Integrator Tstops Tests on the Interval $(["[-1, 0]", "[0, 1]"][i])" for (i, tdir) in enumerate(
        [
            -1.0;
            1.0
        ]
    )
    prob2 = remake(prob_ode_linear, tspan = (0.0, tdir * 1.0))
    integrator = init(prob2, Tsit5())
    tstops = tdir .* [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]
    for tstop in tstops
        add_tstop!(integrator, tstop)
    end
    @test_throws ErrorException add_tstop!(integrator, -0.1 * tdir)
    solve!(integrator)
    for tstop in tstops
        @test tstop ∈ integrator.sol.t
    end
end

@testset "Tstops Eps" begin
    function de(du, u, p, t) # specific DE does not impact the issue
        a, b = p
        du[1] = a * u[1]
        du[2] = b * u[2]
    end

    saveat = [0.0, 0.0094777, 1.5574]
    tstop = 0.010823

    affect!(integrator) = integrator.u[1] += 1.0
    condition(u, t, integrator) = t == tstop
    callback = DiscreteCallback(condition, affect!)

    prob = ODEProblem(de, zeros(2), (-1, 3.0), rand(2))
    sol = solve(prob, Tsit5(), saveat = saveat, tstops = tstop, callback = callback)
    @test sol.t[end] == 1.5574
end

@testset "Tstops Type Conversion" begin
    called = Ref(false)
    tval = rand()
    ff(du, u, p, t) = du .= 0
    cb = DiscreteCallback(
        (u, t, integrator) -> t == Float32(tval),
        integrator -> (called[] = true)
    )
    prob = ODEProblem(ff, [0.0], (0.0f0, 1.0f0))
    sol = solve(prob, Tsit5(), tstops = [tval], callback = cb)
end

@testset "Late binding tstops" begin
    function rhs(u, p, t)
        u * p + t
    end
    prob = ODEProblem(rhs, 1.0, (0.0, 1.0), 0.1; tstops = (p, tspan) -> tspan[1]:p:tspan[2])
    sol = solve(prob, Tsit5())
    @test 0.0:0.1:1.0 ⊆ sol.t
    prob2 = remake(prob; p = 0.07)
    sol2 = solve(prob2, Tsit5())
    @test 0.0:0.07:1.0 ⊆ sol2.t
end

# Tests for issue #2752: tstop overshoot errors with StaticArrays

@testset "StaticArrays vs Arrays with extreme precision" begin
    function precise_dynamics(u, p, t)
        x = @view u[1:2]
        v = @view u[3:4]
        dv = -0.01 * x + 1.0e-6 * sin(100 * t) * SVector{2}(1, 1)
        return SVector{4}(v[1], v[2], dv[1], dv[2])
    end

    function precise_dynamics_array!(du, u, p, t)
        x = @view u[1:2]
        v = @view u[3:4]
        dv = -0.01 * x + 1.0e-6 * sin(100 * t) * [1, 1]
        du[1] = v[1]
        du[2] = v[2]
        du[3] = dv[1]
        du[4] = dv[2]
    end

    u0_static = SVector{4}(1.0, -0.5, 0.01, 0.01)
    u0_array = [1.0, -0.5, 0.01, 0.01]
    tspan = (0.0, 2.0)
    tstops = [0.5, 1.0, 1.5]

    prob_static = ODEProblem(precise_dynamics, u0_static, tspan)
    sol_static = solve(
        prob_static, Vern9(); reltol = 1.0e-12, abstol = 1.0e-15,
        tstops = tstops
    )
    @test SciMLBase.successful_retcode(sol_static)
    for tstop in tstops
        @test tstop ∈ sol_static.t
    end

    prob_array = ODEProblem(precise_dynamics_array!, u0_array, tspan)
    sol_array = solve(
        prob_array, Vern9(); reltol = 1.0e-12, abstol = 1.0e-15,
        tstops = tstops
    )
    @test SciMLBase.successful_retcode(sol_array)
    for tstop in tstops
        @test tstop ∈ sol_array.t
    end

    @test isapprox(sol_static(2.0), sol_array(2.0), rtol = 1.0e-10)
end

@testset "Backward integration with tstop flags" begin
    function decay_ode(u, p, t)
        SA[-0.1 * u[1]]
    end

    u0 = SVector{1}(1.0)
    tspan = (2.0, 0.0)
    tstops = [1.5, 1.0, 0.5]

    prob = ODEProblem(decay_ode, u0, tspan)
    sol = solve(prob, Vern9(); tstops = tstops, reltol = 1.0e-12, abstol = 1.0e-15)
    @test SciMLBase.successful_retcode(sol)
    for tstop in tstops
        @test tstop ∈ sol.t
    end
end

@testset "PresetTimeCallback with tstop flags" begin
    callback_times = Float64[]

    function affect_preset!(integrator)
        push!(callback_times, integrator.t)
        integrator.u += 0.1 * integrator.u
    end

    function simple_growth(u, p, t)
        SA[0.1 * u[1]]
    end

    u0 = SA[1.0]
    tspan = (0.0, 3.0)
    critical_times = [0.5, 1.0, 1.5, 2.0, 2.5]

    preset_cb = PresetTimeCallback(critical_times, affect_preset!)

    prob = ODEProblem(simple_growth, u0, tspan)
    sol = solve(
        prob, Vern9(); tstops = critical_times, callback = preset_cb,
        reltol = 1.0e-10, abstol = 1.0e-12
    )

    @test SciMLBase.successful_retcode(sol)
    for time in critical_times
        @test any(abs.(sol.t .- time) .< 1.0e-10)
    end
    @test length(callback_times) == length(critical_times)
end

@testset "Multiple close tstops with StaticArrays" begin
    function oscillator(u, p, t)
        SVector{2}(u[2], -u[1])
    end

    u0 = SVector{2}(1.0, 0.0)
    tspan = (0.0, 4.0)
    close_tstops = [
        1.0, 1.0 + 1.0e-14, 1.0 + 2.0e-14, 1.0 + 5.0e-14,
        2.0, 2.0 + 1.0e-15, 2.0 + 1.0e-14,
        3.0, 3.0 + 1.0e-13,
    ]

    prob = ODEProblem(oscillator, u0, tspan)
    sol = solve(prob, Vern9(); tstops = close_tstops, reltol = 1.0e-12, abstol = 1.0e-15)
    @test SciMLBase.successful_retcode(sol)
    for time in [1.0, 2.0, 3.0]
        @test any(abs.(sol.t .- time) .< 1.0e-10)
    end
end

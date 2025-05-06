using OrdinaryDiffEq, Test, Random
import ODEProblemLibrary: prob_ode_linear
Random.seed!(100)

@testset "Tstops Tests on the Interval [0, 1]" begin
    prob = prob_ode_linear

    sol = solve(prob, Tsit5(), dt = 1 // 2^(6), tstops = [1 / 2])
    @test 1 // 2 ∈ sol.t

    sol = solve(prob, RK4(), dt = 1 // 3, tstops = [1 / 2], adaptive = false)
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    sol = solve(prob, RK4(), dt = 1 // 3, tstops = [1 / 2],
        d_discontinuities = [-1 / 2, 1 / 2, 3 / 2], adaptive = false)
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    # TODO
    integrator = init(prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4],
        adaptive = false)

    sol = solve(prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4], adaptive = false)
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(prob, RK4(), tstops = [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1],
        adaptive = false)
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(prob, RK4(), tstops = 0:(1 // 16):1, adaptive = false)
    @test sol.t == collect(0:(1 // 16):1)

    sol = solve(prob, RK4(), tstops = range(0, stop = 1, length = 100), adaptive = false)
    @test sol.t == collect(range(0, stop = 1, length = 100))
end

@testset "Integrator Tstops Tests on the Interval $(["[-1, 0]", "[0, 1]"][i])" for (i, tdir) in enumerate([-1.0;
                                                                                                           1.0])
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
    cb = DiscreteCallback((u, t, integrator) -> t == Float32(tval),
        integrator -> (called[] = true))
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

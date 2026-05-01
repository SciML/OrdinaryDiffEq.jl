using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqTsit5
using SciMLBase: ReturnCode
using Test

const prob = prob_dde_constant_2delays_ip

# total order
@testset "total order" begin
    a = Discontinuity(1, 3)

    @test a > 0
    @test a < 2
    @test a == 1

    b = Discontinuity(2.0, 2)

    @test !(a > b)
    @test a < b
    @test a != b

    c = Discontinuity(1.0, 2)

    @test a > c
    @test !(a < c)
    @test a != c
end

# simple DDE example
@testset "DDE" begin
    integrator = init(prob, MethodOfSteps(BS3()))

    # initial discontinuities
    @testset "initial" begin
        @test integrator.tracked_discontinuities == [Discontinuity(0.0, 0)]
        @test length(integrator.d_discontinuities_propagated) == 2 &&
            issubset(
            [Discontinuity(1 / 5, 1), Discontinuity(1 / 3, 1)],
            integrator.d_discontinuities_propagated.valtree
        )
        @test isempty(integrator.opts.d_discontinuities)
        @test isempty(integrator.opts.d_discontinuities_cache)
    end

    # tracked discontinuities
    @testset "tracked" begin
        solve!(integrator)

        discs = [
            Discontinuity(t, order)
                for (t, order) in (
                    (0.0, 0), (1 / 5, 1), (1 / 3, 1), (2 / 5, 2),
                    (8 / 15, 2), (3 / 5, 3),
                    (2 / 3, 2), (11 / 15, 3), (13 / 15, 3),
                )
        ]

        for (tracked, disc) in zip(integrator.tracked_discontinuities, discs)
            @test tracked.t ≈ disc.t && tracked.order == disc.order
        end
    end
end

# additional discontinuities
@testset "DDE with discontinuities" begin
    integrator = init(
        prob, MethodOfSteps(BS3());
        d_discontinuities = [Discontinuity(0.3, 4), Discontinuity(0.6, 5)]
    )

    @test integrator.tracked_discontinuities == [Discontinuity(0.0, 0)]
    @test length(integrator.d_discontinuities_propagated) == 2 &&
        issubset(
        [Discontinuity(1 / 5, 1), Discontinuity(1 / 3, 1)],
        integrator.d_discontinuities_propagated.valtree
    )
    @test length(integrator.opts.d_discontinuities) == 2 &&
        issubset(
        [Discontinuity(0.3, 4), Discontinuity(0.6, 5)],
        integrator.opts.d_discontinuities.valtree
    )
    @test integrator.opts.d_discontinuities_cache ==
        [Discontinuity(0.3, 4), Discontinuity(0.6, 5)]
end

# discontinuities induced by callbacks
@testset "#190" begin
    cb = ContinuousCallback((u, t, integrator) -> 1, integrator -> nothing)
    integrator = init(
        prob_dde_constant_1delay_ip, MethodOfSteps(Tsit5());
        callback = cb
    )
    sol = solve!(integrator)
    for t in 0:5
        @test t ∈ sol.t
        @test Discontinuity(Float64(t), t) ∈ integrator.tracked_discontinuities
    end
end

# OrdinaryDiffEqCore's `shift_past_discontinuity!` (PR #3392) nudges
# `integrator.t` one ULP past each discontinuity. DDE propagates a discontinuity
# at every constant-lag step; without the DDE-side override and the relaxed
# DDE/ODE sync check, the next step throws "cannot advance ODE integrator" as
# soon as integration crosses `t = n * tau`. Regression test for
# SciML/OrdinaryDiffEq.jl#3426.
@testset "constant-lag shift past discontinuity (#3426)" begin
    function bc_model(du, u, h, p, t)
        tau = 1
        d = h(p, t - tau)[3]^2
        du[1] = (1 / (1 + d)) * (0.2 - 0.3) * u[1] - 5 * u[1]
        du[2] = (1 / (1 + d)) * (1 - 0.2 + 0.3) * u[1] +
            (1 / (1 + d)) * (0.2 - 0.3) * u[2] - u[2]
        du[3] = (1 / (1 + d)) * (1 - 0.2 + 0.3) * u[2] - u[3]
        return nothing
    end
    h_3426(p, t) = ones(3)
    prob_3426 = DDEProblem(
        bc_model, [1.0, 1.0, 1.0], h_3426, (0.0, 10.0); constant_lags = [1]
    )
    for alg in (BS3(), Tsit5())
        sol = solve(prob_3426, MethodOfSteps(alg); reltol = 1.0e-7, abstol = 1.0e-10)
        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 10.0
    end
end

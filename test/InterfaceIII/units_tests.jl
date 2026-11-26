using OrdinaryDiffEq, RecursiveArrayTools, Unitful
using LinearAlgebra, Test, ADTypes
using OrdinaryDiffEqExplicitRK, OrdinaryDiffEqFeagin, OrdinaryDiffEqHighOrderRK,
    OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSSPRK
using OrdinaryDiffEqRKN: SharpFineRKN6

@testset "SharpFineRKN6 unitful acceleration" begin
    function damped!(ddu, du, u, p, t)
        @. ddu = -u / (1u"s")^2 - 0.5 * du / (1u"s")
    end
    prob = SecondOrderODEProblem(damped!, [0.0u"m/s"], [1.0u"m"], (0.0u"s", 1.0u"s"))
    sol = solve(prob, SharpFineRKN6(); adaptive = false, dt = 0.125u"s")
    theta = sqrt(15) / 4
    @test successful_retcode(sol)
    @test sol.u[end].x[1][1] ≈ (-4 * exp(-1 / 4) * sin(theta) / sqrt(15))u"m/s" rtol = 1.0e-6
    @test sol.u[end].x[2][1] ≈ (exp(-1 / 4) * (cos(theta) + sin(theta) / sqrt(15)))u"m" rtol = 1.0e-6
end

@testset "Algorithms" begin
    algs = [
        Euler(), Midpoint(), Heun(), Ralston(), RK4(), SSPRK104(), SSPRK22(), SSPRK33(),
        SSPRK43(), SSPRK432(), BS3(), BS5(), DP5(), DP8(), Feagin10(), Feagin12(),
        Feagin14(), TanYam7(), Tsit5(), TsitPap8(), Vern6(), Vern7(), Vern8(), Vern9(),
    ]

    @testset "Scalar units" begin
        f(y, p, t) = 0.5 * y / 3.0u"s"
        u0 = 1.0u"N"
        prob = ODEProblem(f, u0, (0.0u"s", 1.0u"s"))

        @testset "$(typeof(alg))" for alg in algs
            sol = solve(prob, alg, dt = 1u"s" / 10)
            @test successful_retcode(sol)
        end

        sol = solve(prob, ExplicitRK(), dt = 1u"s" / 10)
    end

    @testset "2D units" begin
        f(dy, y, p, t) = (dy .= 0.5 .* y ./ 3.0u"s")
        u0 = [
            1.0u"N" 2.0u"N"
            3.0u"N" 1.0u"N"
        ]
        prob = ODEProblem(f, u0, (0.0u"s", 1.0u"s"))

        @testset "$(typeof(alg))" for alg in algs
            sol = solve(prob, alg, dt = 1u"s" / 10)
            @test successful_retcode(sol)
        end

        sol = solve(prob, ExplicitRK(), dt = 1u"s" / 10)
    end
end

@testset "Mixed units" begin
    @testset "With ArrayPartition" begin
        r0 = [1131.34, -2282.343, 6672.423]u"km"
        v0 = [-5.64305, 4.30333, 2.42879]u"km/s"
        Δt = 86400.0 * 365u"s"
        μ = 398600.4418u"km^3/s^2"
        rv0 = ArrayPartition(r0, v0)

        function f(dy, y, μ, t)
            r = norm(y.x[1])
            dy.x[1] .= y.x[2]
            dy.x[2] .= -μ .* y.x[1] / r^3
        end

        prob = ODEProblem(f, rv0, (0.0u"s", 1.0u"s"), μ)
        @testset "$(typeof(alg))" for alg in [Tsit5()]
            sol = solve(prob, alg)
        end

        @testset "$(typeof(alg))" for alg in [
                AutoVern6(Rodas5P(autodiff = AutoFiniteDiff())),
                AutoVern7(Rodas5P(autodiff = AutoFiniteDiff())),
                AutoVern8(Rodas5P(autodiff = AutoFiniteDiff())),
                AutoVern9(Rodas5P(autodiff = AutoFiniteDiff())),
            ]
            @test_broken sol = solve(prob, alg)
        end
    end

    @testset "Without ArrayPartition" begin
        # coordinate: u = [position, momentum]
        # parameters: p = [mass, force constant]
        function f_harmonic!(du, u, p, t)
            du[1] = u[2] / p[1]
            du[2] = -p[2] * u[1]
        end

        mass = 1.0u"kg"
        k = 1.0u"N/m"
        p = [mass, k]

        u0 = [1.0u"m", 0.0u"kg*m/s"] # initial values (position, momentum)
        tspan = (0.0u"s", 10.0u"s")
        prob = ODEProblem(f_harmonic!, u0, tspan, p)
        @test successful_retcode(solve(prob, Tsit5()))
    end

    @testset "Unitful time with unitless state" begin
        using OrdinaryDiffEq, Unitful
        u0 = 30.0
        tspan = (0.0u"s", 10.0u"s")
        prob1 = ODEProblem((du, u, t, p) -> (du[1] = -0.2u"1/s" * u[1]), [u0], tspan)
        prob2 = ODEProblem((u, t, p) -> (-0.2u"1/s" * u[1]), u0, tspan)
        prob3 = ODEProblem((u, t, p) -> [-0.2u"1/s" * u[1]], [u0], tspan)
        for prob in [prob1, prob2, prob3]
            @test successful_retcode(solve(prob, Tsit5()))
        end
    end
end

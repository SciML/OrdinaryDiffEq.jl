using OrdinaryDiffEq, DelayDiffEq, BoundaryValueDiffEq, DiffEqDevTools, DiffEqBase, Test
using ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear
using DDEProblemLibrary: prob_dde_constant_1delay_ip
using BVProblemLibrary: prob_bvp_linear_1
using Test

## Setup Tests

prob = prob_ode_linear

probs = Vector{ODEProblem}(undef, 2)
probs[1] = prob
probs[2] = prob_ode_2Dlinear

tspan = [0, 1]
tspans = Vector{Vector{Int64}}(undef, 2)
tspans[1] = tspan;
tspans[2] = tspan;

setups = [Dict(:alg => RK4()); Dict(:alg => Euler()); Dict(:alg => BS3());
          Dict(:alg => Midpoint()); Dict(:alg => BS5()); Dict(:alg => DP5())]

t1 = @elapsed sol = solve(prob, RK4(), dt = 1 / 2^(4))
t2 = @elapsed sol2 = solve(prob, setups[1][:alg], dt = 1 / 2^(4))

@test (sol2.u[end] == sol.u[end])

test_sol_2Dlinear = TestSolution(
    solve(prob_ode_2Dlinear, Vern7(), abstol = 1 / 10^14, reltol = 1 / 10^14))

prob_lotka = begin
    function lotka(du, u, p, t)
        du[1] = 1.5 * u[1] - u[1] * u[2]
        du[2] = -3 * u[2] + u[1] * u[2]
    end
    ODEProblem(lotka, [1.0; 1.0], (0.0, 10.0))
end
test_sol_lotka = TestSolution(
    solve(prob_lotka, Vern7(), abstol = 1 / 10^14, reltol = 1 / 10^14))

## Shootout Tests
@testset "Shootout Tests" begin
    println("Shootout Tests")

    shoot = Shootout(prob, setups, dt = 1 / 2^(4))

    #show(shoot)
    #println(shoot)
    shoot[end]
    @test shoot.names == ["RK4", "Euler", "BS3", "Midpoint", "BS5", "DP5"]

    set = ShootoutSet(probs, setups; dt = 1 / 2^(4))

    #println(set[1])
    #println(set[:])
    set[end]
    set[1][:]
    @test all(
        x -> x.names == ["RK4", "Euler", "BS3", "Midpoint", "BS5", "DP5"], set.shootouts)

    @testset "RadauIIA5 and RosShamp4" begin
        __f(u, p, t) = 1.01 * u
        u0 = 1 / 2
        tspan = (0.0, 1.0)
        prob = ODEProblem(__f, u0, tspan)
        sol = solve(prob, Rodas4(), reltol = 1e-8, abstol = 1e-8)
        setups = [Dict(:alg => RadauIIA5()),
            Dict(:alg => RosShamp4())]
        shoot = Shootout(prob, setups; appxsol = TestSolution(sol))
        @test shoot.names == ["RadauIIA5", "RosShamp4"]
    end
end

## WorkPrecision Tests
@testset "WorkPrecision Tests" begin
    println("WorkPrecision Tests")
    @testset "Test DP5 and Tsit5" begin
        println("Test DP5")
        abstols = 1 ./ 10 .^ (3:10)
        reltols = 1 ./ 10 .^ (3:10)
        wp = WorkPrecision(
            prob_ode_linear, DP5(), abstols, reltols; name = "Dormand-Prince 4/5")

        wp[1]
        wp[:]
        wp[end]
        #show(wp)
    end

    @testset "Test setups" begin
        abstols = 1 ./ 10 .^ (3:10)
        reltols = 1 ./ 10 .^ (3:10)
        wp_set = WorkPrecisionSet(
            prob_ode_linear, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)

        wp_set[1]
        wp_set[:]
        wp_set[end]
        #println(wp_set)
        #show(wp_set)
        @test (minimum(diff(wp_set[2].errors.final) .== 0)) # The errors for a fixed timestep method should be constant
        @test wp_set.names == ["RK4", "Euler", "BS3", "Midpoint", "BS5", "DP5"]
    end

    @testset "Test DP5 and Tsit5" begin
        prob = prob_ode_2Dlinear

        abstols = 1 ./ 10 .^ (3:7)
        reltols = 1 ./ 10 .^ (0:4)

        setups = [Dict(:alg => DP5())
                  Dict(:alg => Tsit5(),
                      :abstols => 1 ./ 10 .^ (4:7),
                      :reltols => 1 ./ 10 .^ (1:4))]

        println("Test DP5 and Tsit5")
        wp = WorkPrecisionSet(prob, abstols, reltols, setups; save_everystep = false)

        @test wp.names == ["DP5", "Tsit5"]
        @test length(wp[1]) == 5
        @test length(wp[2]) == 4
    end

    @testset "Test DP5, Tsit5, and Vern6" begin
        abstols = 1 ./ 10 .^ (6:9)
        reltols = 1 ./ 10 .^ (3:6)

        setups = [Dict(:alg => DP5())
                  Dict(:alg => Tsit5())
                  Dict(:alg => Vern6())]
        println("Test DP5, Tsit5, and Vern6")
        wp = WorkPrecisionSet(
            prob_lotka, abstols, reltols, setups; appxsol = test_sol_lotka,
            save_everystep = false, numruns = 20, maxiters = 10000)
        @test wp.names == ["DP5", "Tsit5", "Vern6"]
    end

    # Dual Problem
    @testset "Dual Problem" begin
        probs = [prob_lotka, prob_ode_2Dlinear]
        setups = [Dict(:alg => DP5())
                  Dict(:alg => Tsit5(), :prob_choice => 2)
                  Dict(:alg => Vern6(), :prob_choice => 2)]
        abstols = 1 ./ 10 .^ (6:9)
        reltols = 1 ./ 10 .^ (3:6)
        println("Test DP5, Tsit5, and Vern6")
        wp = WorkPrecisionSet(
            probs, abstols, reltols, setups; appxsol = [test_sol_lotka, nothing],
            save_everystep = false, numruns = 20, maxiters = 10000)
        @test wp.names == ["DP5", "Tsit5", "Vern6"]
        wp = WorkPrecisionSet(
            probs, abstols, reltols, setups; appxsol = [test_sol_lotka, test_sol_2Dlinear],
            save_everystep = false, numruns = 20, maxiters = 10000)
        @test wp.names == ["DP5", "Tsit5", "Vern6"]
    end

    # Dual DAE Problems
    @testset "Dual DAE Problems" begin
        function rober(du, u, p, t)
            y₁, y₂, y₃ = u
            k₁, k₂, k₃ = p
            du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
            du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
            du[3] = k₂ * y₂^2
            nothing
        end
        prob1 = ODEProblem(rober, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])

        ode_ref_sol = solve(prob1, Rodas5(), abstol = 1 / 10^14, reltol = 1 / 10^14)

        function dae_rober(out, du, u, p, t)
            out[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
            out[2] = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
            out[3] = u[1] + u[2] + u[3] - 1.0
        end
        u₀ = [1.0, 0, 0]
        du₀ = [-0.04, 0.04, 0.0]
        tspan = (0.0, 100000.0)

        differential_vars = [true, true, false]
        prob2 = DAEProblem(dae_rober, du₀, u₀, tspan, differential_vars = differential_vars)

        dae_ref_sol = solve(prob2, DFBDF(), abstol = 1 / 10^14, reltol = 1 / 10^14)

        probs = [prob1, prob2]
        setups = [Dict(:alg => Rodas5())
                  Dict(:alg => DFBDF(), :prob_choice => 2)]
        abstols = 1 ./ 10 .^ (6:9)
        reltols = 1 ./ 10 .^ (3:6)
        wp = WorkPrecisionSet(
            probs, abstols, reltols, setups; appxsol = [ode_ref_sol, dae_ref_sol],
            save_everystep = false, numruns = 20, maxiters = 10000)
        @test wp.names == ["Rodas5", "DFBDF"]
    end

    # DDE problem
    @testset "DDE problem" begin
        prob = prob_dde_constant_1delay_ip

        abstols = 1 ./ 10 .^ (7:10)
        reltols = 1 ./ 10 .^ (4:7)
        sol = solve(
            prob, MethodOfSteps(Vern9(), fpsolve = NLFunctional(; max_iter = 1000));
            reltol = 1e-8, abstol = 1e-8)
        test_sol = TestSolution(sol)

        setups = [Dict(:alg => MethodOfSteps(BS3()))
                  Dict(:alg => MethodOfSteps(Tsit5()))]
        println("Test MethodOfSteps BS3 and Tsit5")
        wp = WorkPrecisionSet(prob, abstols, reltols, setups; appxsol = test_sol)
        @test wp.names == ["BS3", "Tsit5"]
        println("DDE Done")
    end

    # BVP problem
    @testset "BVP Problem" begin
        prob = prob_bvp_linear_1
        abstols = 1.0 ./ 10.0 .^ (2:3)
        reltols = 1.0 ./ 10.0 .^ (2:3)
        sol = solve(prob, Shooting(Tsit5()), abstol = 1e-14, reltol = 1e-14)
        test_sol = TestSolution(sol.t, sol.u)

        setups = [Dict(:alg => MIRK4(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))
                  Dict(:alg => MIRK5(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1))]
        labels = ["MIRK4", "MIRK5"]

        println("Test MIRK4 and MIRK5")
        wp = WorkPrecisionSet(prob,
            abstols,
            reltols,
            setups;
            appxsol = test_sol,
            names = labels,
            print_names = true)
        @test wp.names == ["MIRK4", "MIRK5"]
        println("BVP Done")
    end
end

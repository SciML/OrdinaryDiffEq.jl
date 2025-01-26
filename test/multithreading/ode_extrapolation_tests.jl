# Import packages
using OrdinaryDiffEqExtrapolation, DiffEqDevTools, Test, Random

println("Running on $(Threads.nthreads()) thread(s).")

# Define test problems
# Note that the time span in ODEProblemLibrary is given by
# Float64 numbers

linear = (u, p, t) -> (p * u)
linear_analytic = (u0, p, t) -> u0 * exp(p * t)
prob_ode_bigfloatlinear = ODEProblem(ODEFunction(linear, analytic = linear_analytic),
    big"0.5", (big"0.0", big"1.0"), big"1.01")

f_2dlinear = (du, u, p, t) -> (@. du = p * u)
f_2dlinear_analytic = (u0, p, t) -> @. u0 * exp(p * t)
prob_ode_bigfloat2Dlinear = ODEProblem(
    ODEFunction(f_2dlinear,
        analytic = f_2dlinear_analytic),
    rand(BigFloat, (4, 2)), (big"0.0", big"1.0"),
    big"1.01")

# Prepare tests
Random.seed!(100)
problem_array = [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
dts = 1 .// 2 .^ (8:-1:1)

testTol = 0.2

@testset "Testing extrapolation methods" begin

    # Test AitkenNeville
    @testset "Testing AitkenNeville" begin
        @testset "Testing sequential AitkenNeville" begin
            for prob in problem_array
                global dts

                #  Convergence test
                for j in 1:4
                    sim = test_convergence(dts, prob,
                        AitkenNeville(max_order = j,
                            min_order = j, init_order = j,
                            threading = false))
                    @test sim.𝒪est[:final]≈j atol=testTol
                end

                # Regression test
                sol = solve(prob,
                    AitkenNeville(max_order = 9, min_order = 1,
                        init_order = 9, threading = false), reltol = 1e-3)
                @test length(sol.u) < 15
                sol = solve(prob,
                    AitkenNeville(max_order = 9, min_order = 1,
                        init_order = 9, threading = false), reltol = 1e-6)
                @test length(sol.u) < 18
            end
        end
        @testset "Testing threaded AitkenNeville" begin
            for prob in problem_array
                global dts

                #  Convergence test
                for j in 1:4
                    sim = test_convergence(dts, prob,
                        AitkenNeville(max_order = j,
                            min_order = j, init_order = j,
                            threading = true))
                    @test sim.𝒪est[:final]≈j atol=testTol
                end

                # Regression test
                sol = solve(prob,
                    AitkenNeville(max_order = 9, min_order = 1,
                        init_order = 9, threading = true), reltol = 1e-3)
                @test length(sol.u) < 15
                sol = solve(prob,
                    AitkenNeville(max_order = 9, min_order = 1,
                        init_order = 9, threading = true), reltol = 1e-6)
                @test length(sol.u) < 18
            end
        end
    end # AitkenNeville

    # Define the subdividing sequences
    sequence_array = [:harmonic, :romberg, :bulirsch]

    @testset "Testing ImplicitEulerExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            newTol = 0.35
            #  Convergence test
            for j in 1:4
                alg = ImplicitEulerExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq, threading = false)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈alg.init_order + 1.1 atol=newTol #Superconvergence
            end
            # Regression test
            sol = solve(prob,
                ImplicitEulerExtrapolation(max_order = 9, min_order = 1,
                    init_order = 9, sequence = seq,
                    threading = false), reltol = 1e-3)
            @test length(sol.u) < 15
        end
    end

    @testset "Testing ImplicitEulerExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            newTol = 0.35
            #  Convergence test
            for j in 1:4
                alg = ImplicitEulerExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq, threading = true)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈alg.init_order + 1.1 atol=newTol #Superconvergence
            end
            # Regression test
            sol = solve(prob,
                ImplicitEulerExtrapolation(max_order = 9, min_order = 1,
                    init_order = 9, sequence = seq,
                    threading = true), reltol = 1e-3)
            @test length(sol.u) < 15
        end
    end

    @testset "Testing ImplicitEulerBarycentricExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            newTol = 0.35
            #  Convergence test
            for j in 1:4
                alg = ImplicitEulerBarycentricExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq,
                    threading = false)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈alg.init_order + 0.5 atol=newTol #Superconvergence
            end
            # Regression test
            sol = solve(prob,
                ImplicitEulerBarycentricExtrapolation(max_order = 9, min_order = 1,
                    init_order = 9,
                    sequence = seq,
                    threading = false),
                reltol = 1e-3)
            @test length(sol.u) < 15
        end
    end

    @testset "Testing ImplicitEulerBarycentricExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            newTol = 0.35
            #  Convergence test
            for j in 1:4
                alg = ImplicitEulerBarycentricExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq,
                    threading = true)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈alg.init_order + 0.5 atol=newTol #Superconvergence
                algp = ImplicitEulerBarycentricExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq,
                    threading = OrdinaryDiffEqExtrapolation.PolyesterThreads())
                simp = test_convergence(dts, prob, algp)
                @test simp.𝒪est[:final]≈algp.init_order + 0.5 atol=newTol #Superconvergence
            end
            # Regression test
            sol = solve(prob,
                ImplicitEulerBarycentricExtrapolation(max_order = 9, min_order = 1,
                    init_order = 9,
                    sequence = seq,
                    threading = true),
                reltol = 1e-3)
            @test length(sol.u) < 15
        end
    end

    @testset "Testing ImplicitDeuflhardExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            # Convergence test
            for j in 1:6
                alg = ImplicitDeuflhardExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq,
                    threading = OrdinaryDiffEqExtrapolation.Sequential())
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
            end

            # Regression test
            alg = ImplicitDeuflhardExtrapolation(max_order = 9, min_order = 1,
                init_order = 9, sequence = seq,
                threading = false)
            sol = solve(prob, alg, reltol = 1e-3)
            @test length(sol.u) < 10
        end
    end

    @testset "Testing ImplicitDeuflhardExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            # Convergence test
            for j in 1:6
                alg = ImplicitDeuflhardExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq, threading = true)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
            end

            # Regression test
            alg = ImplicitDeuflhardExtrapolation(max_order = 9, min_order = 1,
                init_order = 9, sequence = seq,
                threading = OrdinaryDiffEqExtrapolation.BaseThreads())
            sol = solve(prob, alg, reltol = 1e-3)
            @test length(sol.u) < 10
        end
    end

    @testset "Testing ImplicitHairerWannerExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            # Convergence test
            for j in 1:6
                alg = ImplicitHairerWannerExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq, threading = false)
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) - 1 atol=testTol
            end

            alg = ImplicitHairerWannerExtrapolation(max_order = 9, min_order = 1,
                init_order = 9, sequence = seq,
                threading = false)
            sol = solve(prob, alg, reltol = 1e-3)
            @test length(sol.u) < 10
        end
    end

    @testset "Testing ImplicitHairerWannerExtrapolation" begin
        for prob in problem_array,
            seq in sequence_array

            global dts

            # Convergence test
            for j in 1:6
                alg = ImplicitHairerWannerExtrapolation(min_order = j,
                    init_order = j, max_order = j,
                    sequence = seq,
                    threading = OrdinaryDiffEqExtrapolation.PolyesterThreads())
                sim = test_convergence(dts, prob, alg)
                @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) - 1 atol=testTol
            end

            alg = ImplicitHairerWannerExtrapolation(max_order = 9, min_order = 1,
                init_order = 9, sequence = seq,
                threading = true)
            sol = solve(prob, alg, reltol = 1e-3)
            @test length(sol.u) < 10
        end
    end

    # Test ExtrapolationMidpointDeuflhard
    @testset "Testing ExtrapolationMidpointDeuflhard" begin
        @testset "Testing sequential ExtrapolationMidpointDeuflhard" begin
            for prob in problem_array,
                seq in sequence_array

                global dts

                # Convergence test
                for j in 1:6
                    alg = ExtrapolationMidpointDeuflhard(min_order = j,
                        init_order = j, max_order = j,
                        sequence = seq,
                        threading = OrdinaryDiffEqExtrapolation.Sequential())
                    sim = test_convergence(dts, prob, alg)
                    @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
                end

                # Regression test
                alg = ExtrapolationMidpointDeuflhard(max_order = 9, min_order = 1,
                    init_order = 9, sequence = seq,
                    threading = false)
                sol = solve(prob, alg, reltol = 1e-3)
                @test length(sol.u) < 10
            end
        end
        @testset "Testing threaded ExtrapolationMidpointDeuflhard" begin
            for prob in problem_array,
                seq in sequence_array

                global dts

                # Convergence test
                for j in 1:6
                    alg = ExtrapolationMidpointDeuflhard(min_order = j,
                        init_order = j, max_order = j,
                        sequence = seq, threading = true)
                    sim = test_convergence(dts, prob, alg)
                    @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
                end

                # Regression test
                alg = ExtrapolationMidpointDeuflhard(max_order = 9, min_order = 1,
                    init_order = 9, sequence = seq,
                    threading = true)
                sol = solve(prob, alg, reltol = 1e-3)
                @test length(sol.u) < 10
            end
        end
    end # ExtrapolationMidpointDeuflhard

    # Test ExtrapolationMidpointHairerWanner
    @testset "Testing ExtrapolationMidpointHairerWanner" begin
        @testset "Testing sequential ExtrapolationMidpointHairerWanner" begin
            for prob in problem_array,
                seq in sequence_array

                global dts

                # Convergence test
                for j in 1:6
                    alg = ExtrapolationMidpointHairerWanner(min_order = j,
                        init_order = j, max_order = j,
                        sequence = seq,
                        threading = false)
                    sim = test_convergence(dts, prob, alg)
                    @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
                end

                # Regression test
                alg = ExtrapolationMidpointHairerWanner(max_order = 9, min_order = 2,
                    init_order = 9, sequence = seq,
                    threading = false)
                sol = solve(prob, alg, reltol = 1e-3)
                @test length(sol.u) < 10
            end
        end
        @testset "Testing threaded ExtrapolationMidpointHairerWanner" begin
            for prob in problem_array,
                seq in sequence_array

                global dts

                # Convergence test
                for j in 1:6
                    alg = ExtrapolationMidpointHairerWanner(min_order = j,
                        init_order = j, max_order = j,
                        sequence = seq,
                        threading = true)
                    sim = test_convergence(dts, prob, alg)
                    @test sim.𝒪est[:final]≈2 * (alg.init_order + 1) atol=testTol
                end

                # Regression test
                alg = ExtrapolationMidpointHairerWanner(max_order = 9, min_order = 2,
                    init_order = 9, sequence = seq,
                    threading = true)
                sol = solve(prob, alg, reltol = 1e-3)
                @test length(sol.u) < 10
            end
        end
    end # ExtrapolationMidpointHairerWanner

    @testset "Regression Test Float32 and Float64 Fallbacks" begin
        prob_ode_2Dlinear = ODEProblem(
            ODEFunction(f_2dlinear,
                analytic = f_2dlinear_analytic),
            Float64.(prob_ode_bigfloat2Dlinear.u0), (0.0, 1.0),
            1.01)
        s1 = solve(prob_ode_bigfloat2Dlinear, ExtrapolationMidpointDeuflhard())
        s2 = solve(prob_ode_2Dlinear, ExtrapolationMidpointDeuflhard())
        @test all(all(s1[i] - s2[i] .< 5e-14) for i in 1:length(s1))

        prob_ode_2Dlinear = ODEProblem(
            ODEFunction(f_2dlinear,
                analytic = f_2dlinear_analytic),
            Float32.(prob_ode_bigfloat2Dlinear.u0),
            (0.0f0, 1.0f0), 1.01f0)
        s1 = solve(prob_ode_bigfloat2Dlinear, ExtrapolationMidpointDeuflhard())
        s2 = solve(prob_ode_2Dlinear, ExtrapolationMidpointDeuflhard())
        @test all(all(s1[i] - s2[i] .< 5e-6) for i in 1:length(s1))
    end
end # Extrapolation methods

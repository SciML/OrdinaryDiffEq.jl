using OrdinaryDiffEqStabilizedRK, DiffEqDevTools, Test, LinearAlgebra, Random
using OrdinaryDiffEqStabilizedRK: maxeig!
using SciMLBase: ReturnCode, step!
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
probArr = Vector{ODEProblem}(undef, 2)
probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
    Random.seed!(123)
    eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
    for iip in [true, false], alg in [ROCK4(), ROCK4(eigen_est = eigen_est)]
        println(typeof(alg))
        A = randn(20, 20)
        test_f(u, p, t) = A * u
        test_f(du, u, p, t) = mul!(du, A, u)
        prob = ODEProblem{iip}(test_f, randn(20), (0, 1.0))
        integrator = init(prob, alg)
        eigm = maximum(abs.(eigvals(A)))
        maxeig!(integrator, integrator.cache)
        eigest = integrator.eigen_est
        @test eigest ≈ eigm rtol = 0.35

        A = A - 1.0e2I
        test_stiff(u, p, t) = A * u
        test_stiff(du, u, p, t) = mul!(du, A, u)
        prob = ODEProblem{iip}(test_stiff, ones(20), (0, 1.0))
        @test_nowarn solve(prob, alg)
    end

    # RKL methods use the same path for the eigenvalue tieratin
    @testset "RKL power iteration" begin
        Random.seed!(456)
        eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
        for iip in [true, false], alg in [
                    RKL2(), RKL2(eigen_est = eigen_est),
                    RKL1(), RKL1(eigen_est = eigen_est),
                ]
            println(typeof(alg))
            A = randn(20, 20)
            A = A - maximum(abs.(eigvals(A))) * I   # negative definite to test parabolic stability
            test_f(u, p, t) = A * u
            test_f(du, u, p, t) = mul!(du, A, u)
            prob = ODEProblem{iip}(test_f, randn(20), (0.0, 1.0))
            integrator = init(prob, alg)
            eigm = maximum(abs.(eigvals(A)))
            maxeig!(integrator, integrator.cache)
            eigest = integrator.eigen_est
            @test eigest ≈ eigm rtol = 0.35
        end
    end
end

@testset "Runge-Kutta-Chebyshev and Runge-Kutta-Legendre Convergence Tests" begin
    dts = 1 .// 2 .^ (8:-1:4)
    testTol = 0.1
    for prob in probArr
        println("ROCK2")
        #default ROCK2
        sim = test_convergence(dts, prob, ROCK2())
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        #testing ROCK2 for different minimum stages to insure that the constants are right
        sim = test_convergence(dts, prob, ROCK2(min_stages = 5))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, ROCK2(min_stages = 10))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, ROCK2(min_stages = 21))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        #default ROCK4
        println("ROCK4")
        sim = test_convergence(dts, prob, ROCK4())
        @test sim.𝒪est[:l∞] ≈ 4 atol = testTol
        #testing ROCK4 for different minimum stages to insure that the constants are right
        sim = test_convergence(dts, prob, ROCK4(min_stages = 6))
        @test sim.𝒪est[:l∞] ≈ 4 atol = testTol
        sim = test_convergence(dts, prob, ROCK4(min_stages = 10))
        @test sim.𝒪est[:l∞] ≈ 4 atol = testTol
        sim = test_convergence(dts, prob, ROCK4(min_stages = 21))
        @test sim.𝒪est[:l∞] ≈ 4 atol = testTol

        println("RKC")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKC(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKC(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 10000 / integrator.dt
        sim = test_convergence(dts, prob, RKC(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        println("TSRKC3")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC3(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 3 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC3(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 3 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 10000 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC3(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 3 atol = testTol
        println("SERK2")
        sim = test_convergence(dts, prob, SERK2())
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        println("ESERK4")
        sim = test_convergence(dts, prob, ESERK4())
        @test sim.𝒪est[:l∞] ≈ 4 atol = testTol

        println("RKL1")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKL1(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKL1(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

        println("RKL2")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKL2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKL2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKL2(min_stages = 3, eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, RKL2(min_stages = 5, eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, RKL2(min_stages = 11, eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
    end
    dts = 1 .// 2 .^ (6:-1:2)
    for prob in probArr
        println("ESERK5")
        sim = test_convergence(dts, prob, ESERK5())
        @test sim.𝒪est[:l∞] ≈ 5 atol = testTol
    end
end

@testset "Number of function evaluations" begin
    x = Ref(0)
    u0 = [1.0, 1.0]
    tspan = (0.0, 1.0)
    probop = ODEProblem(u0, tspan) do u, p, t
        x[] += 1
        return -500 * u
    end
    probip = ODEProblem(u0, tspan) do du, u, p, t
        x[] += 1
        @. du = -500 * u
        return nothing
    end

    @testset "$prob" for prob in [probop, probip]
        eigen_est = (integrator) -> integrator.eigen_est = 500
        algs = [
            ROCK2(), ROCK2(eigen_est = eigen_est),
            ROCK4(), ROCK4(eigen_est = eigen_est),
            RKC(), RKC(eigen_est = eigen_est),
            TSRKC3(), TSRKC3(eigen_est = eigen_est),
            SERK2(), SERK2(eigen_est = eigen_est),
            ESERK4(), ESERK4(eigen_est = eigen_est),
            ESERK5(), ESERK5(eigen_est = eigen_est),
            RKL1(), RKL1(eigen_est = eigen_est),
            RKL2(), RKL2(eigen_est = eigen_est),
        ]
        @testset "$alg" for alg in algs
            x[] = 0
            sol = solve(prob, alg)
            @test x[] == sol.stats.nf
        end
    end
end

@testset "Allocations" begin
    u0 = [1.0, 1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(u0, tspan) do du, u, p, t
        @. du = -500 * u
        return nothing
    end

    eigen_est = (integrator) -> integrator.eigen_est = 500
    algs = [
        ROCK2(), ROCK2(eigen_est = eigen_est),
        ROCK4(), ROCK4(eigen_est = eigen_est),
        RKC(), RKC(eigen_est = eigen_est),
        TSRKC3(), TSRKC3(eigen_est = eigen_est),
        SERK2(), SERK2(eigen_est = eigen_est),
        ESERK4(), ESERK4(eigen_est = eigen_est),
        ESERK5(), ESERK5(eigen_est = eigen_est),
        RKL1(), RKL1(eigen_est = eigen_est),
        RKL2(), RKL2(eigen_est = eigen_est),
    ]
    @testset "$alg" for alg in algs
        # compile once
        integrator = init(prob, alg; save_everystep = false)
        solve!(integrator)
        # check allocations
        integrator = init(prob, alg; save_everystep = false)
        allocs = @allocations solve!(integrator)
        # sizehint! calls in _postamble! add allocations to release excess memory
        expected_allocs = VERSION >= v"1.11" ? 7 : 6
        @test allocs <= expected_allocs
    end
end

# RKL specific tests
@testset "RKL-specific tests" begin

    @testset "Constructor enforces odd stages" begin
        @test RKL1(min_stages = 4).min_stages == 5
        @test RKL2(min_stages = 4).min_stages == 5
        @test RKL1(min_stages = 3).min_stages == 3
        @test RKL2(min_stages = 3).min_stages == 3
        @test RKL1(max_stages = 200).max_stages == 199
        @test RKL2(max_stages = 200).max_stages == 199
        @test RKL1(max_stages = 199).max_stages == 199
        @test RKL2(max_stages = 199).max_stages == 199
        @test RKL1(min_stages = 7, max_stages = 5).max_stages >= 7
        @test RKL2(min_stages = 7, max_stages = 5).max_stages >= 7
    end


    @testset "Stiff Parabolic Problem" begin
        N = 20
        α = 0.1
        dx = 1.0 / (N + 1)
        xs = collect(range(dx, 1 - dx; length = N))
        D2 = Tridiagonal(fill(α / dx^2, N - 1), fill(-2α / dx^2, N), fill(α / dx^2, N - 1))
        u0 = sin.(π .* xs)
        tspan = (0.0, 0.5)
        prob = ODEProblem((du, u, p, t) -> mul!(du, D2, u), u0, tspan)

        for alg in [RKL1(), RKL2()]
            sol = solve(prob, alg, abstol = 1.0e-4, reltol = 1.0e-4)
            @test sol.retcode == ReturnCode.Success
            exact = exp(-π^2 * α * tspan[2]) .* sin.(π .* xs)
            @test norm(sol.u[end] - exact) < 0.02
        end
    end

    @testset "Odd Stage Counts" begin
        stages_seen = Int[]
        u0 = ones(5)
        tspan = (0.0, 0.1)
        prob = ODEProblem((du, u, p, t) -> (@. du = -100 * u), u0, tspan)
        eigen_est = (integrator) -> (integrator.eigen_est = 100.0)
        integrator = init(prob, RKL2(eigen_est = eigen_est); save_everystep = false)
        while integrator.t < tspan[2] - 100 * eps(typeof(tspan[2]))
            step!(integrator)
            cc = integrator.cache
            mdeg = hasproperty(cc, :constantcache) ? cc.constantcache.mdeg : cc.mdeg
            push!(stages_seen, mdeg)
        end
        @test all(isodd, stages_seen)
        @test all(s -> s >= 3, stages_seen)
    end
end

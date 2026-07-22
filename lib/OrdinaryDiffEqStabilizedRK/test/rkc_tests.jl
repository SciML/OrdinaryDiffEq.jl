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

    @testset "RKG power iteration" begin
        Random.seed!(789)
        eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
        for iip in [true, false], alg in [
                    RKG2(), RKG2(eigen_est = eigen_est),
                    RKG1(), RKG1(eigen_est = eigen_est),
                ]
            println(typeof(alg))
            A = randn(20, 20)
            A = A - maximum(abs.(eigvals(A))) * I
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

@testset "Runge-Kutta-Chebyshev, Runge-Kutta-Legendre, and Runge-Kutta-Gegenbauer Convergence Tests" begin
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
        println("RKMC2")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKMC2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKMC2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 10000 / integrator.dt
        sim = test_convergence(dts, prob, RKMC2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        println("TSRKC2")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 10000 / integrator.dt
        sim = test_convergence(dts, prob, TSRKC2(eigen_est = eigen_est))
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

        println("RKG1")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKG1(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKG1(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

        println("RKG2")
        eigen_est = (integrator) -> integrator.eigen_est = 1 / integrator.dt
        sim = test_convergence(dts, prob, RKG2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        eigen_est = (integrator) -> integrator.eigen_est = 100 / integrator.dt
        sim = test_convergence(dts, prob, RKG2(eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, RKG2(min_stages = 3, eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
        sim = test_convergence(dts, prob, RKG2(min_stages = 5, eigen_est = eigen_est))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
    end
    dts = 1 .// 2 .^ (6:-1:2)
    for prob in probArr
        println("ESERK5")
        sim = test_convergence(dts, prob, ESERK5())
        @test sim.𝒪est[:l∞] ≈ 5 atol = testTol
    end
end

@testset "In-place and out-of-place agree in the stiff high-stage regime" begin
    # Stiff forced heat equation stepped at fixed dt with dt*λ ≈ 217 so that the
    # stage counts are high enough to exercise the deep recurrences (e.g. SERK2's
    # inner stage loop only runs once mdeg >= 20).
    N = 32
    dx = 1.0 / (N + 1)
    xs = collect(range(dx, 1 - dx; length = N))
    A = Tridiagonal(fill(1 / dx^2, N - 1), fill(-2 / dx^2, N), fill(1 / dx^2, N - 1))
    lam = 4 / dx^2
    g(t) = cos(10t)
    f_oop(u, p, t) = A * u .+ g(t)
    f_iip(du, u, p, t) = (mul!(du, A, u); du .+= g(t); nothing)
    u0 = sin.(pi .* xs)
    tspan = (0.0, 0.5)
    prob_oop = ODEProblem{false}(f_oop, u0, tspan)
    prob_iip = ODEProblem{true}(f_iip, u0, tspan)
    eigen_est = (integrator) -> integrator.eigen_est = lam

    algs = [
        ROCK2(eigen_est = eigen_est), ROCK4(eigen_est = eigen_est),
        RKC(eigen_est = eigen_est), RKMC2(eigen_est = eigen_est),
        TSRKC2(eigen_est = eigen_est), TSRKC3(eigen_est = eigen_est),
        SERK2(eigen_est = eigen_est), ESERK4(eigen_est = eigen_est),
        ESERK5(eigen_est = eigen_est), RKL1(eigen_est = eigen_est),
        RKL2(eigen_est = eigen_est), RKG1(eigen_est = eigen_est),
        RKG2(eigen_est = eigen_est),
    ]
    @testset "$alg" for alg in algs
        sol_oop = solve(
            prob_oop, alg, adaptive = false, dt = 0.05,
            save_everystep = false
        )
        sol_iip = solve(
            prob_iip, alg, adaptive = false, dt = 0.05,
            save_everystep = false
        )
        @test norm(sol_oop.u[end] .- sol_iip.u[end], Inf) < 1.0e-6
    end
end

@testset "RKL non-autonomous convergence" begin
    # Prothero-Robinson-type problem: uᵢ' = -λᵢ*(uᵢ - cos(t)) - sin(t), exact uᵢ = cos(t).
    # Catches wrong or missing stage times, which autonomous problems cannot see.
    lams = [1.0, 5.0, 20.0, 50.0]
    f_na_oop(u, p, t) = @. -p * (u - cos(t)) - sin(t)
    f_na_iip(du, u, p, t) = (@. du = -p * (u - cos(t)) - sin(t); nothing)
    na_analytic(u0, p, t) = cos(t) .* ones(length(p))
    prob_na_oop = ODEProblem(
        ODEFunction{false}(f_na_oop, analytic = na_analytic), ones(4), (0.0, 1.0), lams
    )
    prob_na_iip = ODEProblem(
        ODEFunction{true}(f_na_iip, analytic = na_analytic), ones(4), (0.0, 1.0), lams
    )
    # pure quadrature: u' = cos(t), any correct second-order method gets order 2
    fq_oop(u, p, t) = fill(cos(t), 2)
    fq_iip(du, u, p, t) = (du .= cos(t); nothing)
    q_analytic(u0, p, t) = u0 .+ sin(t)
    prob_q_oop = ODEProblem(
        ODEFunction{false}(fq_oop, analytic = q_analytic), zeros(2), (0.0, 1.0)
    )
    prob_q_iip = ODEProblem(
        ODEFunction{true}(fq_iip, analytic = q_analytic), zeros(2), (0.0, 1.0)
    )

    dts = 1 .// 2 .^ (10:-1:5)
    testTol = 0.25
    # fixed estimate exercises the classical dt -> 0 regime; the 1/dt-scaled
    # estimate keeps the stage count constant (stiff regime)
    fixed_eig = (integrator) -> integrator.eigen_est = 55.0
    scaled_eig = (integrator) -> integrator.eigen_est = 500 / abs(integrator.dt)
    for prob in [prob_na_oop, prob_na_iip, prob_q_oop, prob_q_iip],
            eig in [fixed_eig, scaled_eig]

        sim = test_convergence(dts, prob, RKL1(eigen_est = eig))
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol
        sim = test_convergence(dts, prob, RKL2(eigen_est = eig))
        @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
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
            RKMC2(), RKMC2(eigen_est = eigen_est),
            TSRKC2(), TSRKC2(eigen_est = eigen_est),
            TSRKC3(), TSRKC3(eigen_est = eigen_est),
            SERK2(), SERK2(eigen_est = eigen_est),
            ESERK4(), ESERK4(eigen_est = eigen_est),
            ESERK5(), ESERK5(eigen_est = eigen_est),
            RKL1(), RKL1(eigen_est = eigen_est),
            RKL2(), RKL2(eigen_est = eigen_est),
            RKG1(), RKG1(eigen_est = eigen_est),
            RKG2(), RKG2(eigen_est = eigen_est),
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
        RKMC2(), RKMC2(eigen_est = eigen_est),
        TSRKC2(), TSRKC2(eigen_est = eigen_est),
        TSRKC3(), TSRKC3(eigen_est = eigen_est),
        SERK2(), SERK2(eigen_est = eigen_est),
        ESERK4(), ESERK4(eigen_est = eigen_est),
        ESERK5(), ESERK5(eigen_est = eigen_est),
        RKL1(), RKL1(eigen_est = eigen_est),
        RKL2(), RKL2(eigen_est = eigen_est),
        RKG1(), RKG1(eigen_est = eigen_est),
        RKG2(), RKG2(eigen_est = eigen_est),
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

@testset "RKG-specific tests" begin
    @testset "RKG on stiff parabolic problem" begin
        N = 20
        α = 0.1
        dx = 1.0 / (N + 1)
        xs = collect(range(dx, 1 - dx; length = N))
        D2 = Tridiagonal(
            fill(α / dx^2, N - 1),
            fill(-2α / dx^2, N),
            fill(α / dx^2, N - 1)
        )
        u0 = sin.(π .* xs)
        tspan = (0.0, 0.5)
        prob = ODEProblem((du, u, p, t) -> mul!(du, D2, u), u0, tspan)
        exact = exp(-π^2 * α * tspan[2]) .* sin.(π .* xs)
        for alg in [RKG1(), RKG2()]
            sol = solve(prob, alg, abstol = 1.0e-4, reltol = 1.0e-4)
            @test sol.retcode == ReturnCode.Success
            @test norm(sol.u[end] - exact) < 0.05
        end
    end

    @testset "RKG monotone near Dirichlet boundary" begin
        # Looking for bounded solution that is non-negative on a heat diffusion problem with Dirichlet BCs
        N = 50
        α = 1.166
        dx = 1.0 / (N + 1)
        xs = collect(range(dx, 1 - dx; length = N))

        # spike near left boundary, zero elsewhere
        u0 = zeros(N)
        u0[2] = 100.0

        # Dirichlet BCs
        D2 = Tridiagonal(
            fill(α / dx^2, N - 1),
            fill(-2α / dx^2, N),
            fill(α / dx^2, N - 1)
        )

        tspan = (0.0, 1.0e-4)
        prob = ODEProblem((du, u, p, t) -> mul!(du, D2, u), u0, tspan)

        # RKG2 should remain non negative due to Convex Monotone Property even with Dirichlet BCs
        sol_rkg = solve(prob, RKG2(), abstol = 1.0e-8, reltol = 1.0e-8)
        @test sol_rkg.retcode == ReturnCode.Success
        @test all(sol_rkg.u[end] .>= -1.0e-10)   # allow tiny floating point negatives
        sol_rkl = solve(prob, RKL2(), abstol = 1.0e-8, reltol = 1.0e-8)
        @test sol_rkl.retcode == ReturnCode.Success
    end
end

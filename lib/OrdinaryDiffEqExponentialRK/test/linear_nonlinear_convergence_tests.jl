using OrdinaryDiffEqExponentialRK, Test, DiffEqDevTools, Random, LinearAlgebra, LinearSolve
using OrdinaryDiffEqVerner, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: alg_order
using OrdinaryDiffEqExponentialRK: expRK_operators, phi
using SciMLOperators: MatrixOperator, ScalarOperator

@testset "Caching Out-of-place" begin
    println("Caching Out-of-place")
    μ = 1.01
    linnonlin_f2 = (u, p, t) -> μ * u
    linnonlin_f1 = ScalarOperator(μ)
    linnonlin_fun = SplitFunction(
        linnonlin_f1, linnonlin_f2;
        analytic = (u0, p, t) -> u0 .* exp.(2μ * t)
    )
    prob = SplitODEProblem(linnonlin_fun, 1 / 2, (0.0, 1.0))

    Random.seed!(100)
    dts = 1 ./ 2 .^ (7:-1:4) #14->7 good plot
    for Alg in [
            LawsonEuler,
            NorsettEuler,
            ETDRK2,
            ETDRK3,
            ETDRK4,
            HochOst4,
            Friedli,
            ERK32ZB,
            ERK43ZB,
            ETD2,
            KenCarp3,
            CFNLIRK3,
        ]
        sim = test_convergence(dts, prob, Alg())
        @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.2
    end

    # Dense test
    sim = test_convergence(dts, prob, ETDRK4(), dense_errors = true)
    @test abs(sim.𝒪est[:l2] - 4) < 0.2
end

@testset "Caching Inplace" begin
    println("Caching Inplace")
    μ = 1.01
    u0 = rand(2)
    A = [2.0 -1.0; -1.0 2.0]
    linnonlin_f1 = MatrixOperator(A)
    linnonlin_f2 = (du, u, p, t) -> du .= μ .* u
    linnonlin_fun_iip = SplitFunction(
        linnonlin_f1, linnonlin_f2;
        analytic = (u0, p, t) -> exp((A + μ * I) * t) * u0
    )
    prob = SplitODEProblem(linnonlin_fun_iip, u0, (0.0, 1.0))

    dts = 1 ./ 2 .^ (8:-1:4)
    for Alg in [
            LawsonEuler(),
            NorsettEuler(),
            ETDRK2(),
            ETDRK3(),
            ETDRK4(),
            HochOst4(),
            Friedli(),
            ERK32ZB(),
            ERK43ZB(),
            ETD2(),
        ]
        sim = test_convergence(dts, prob, Alg)
        @test sim.𝒪est[:l2] ≈ alg_order(Alg) atol = 0.15
    end

    dts = 1 ./ 2 .^ (13:-1:9)
    Alg = KenCarp3(linsolve = LinearSolve.KrylovJL_GMRES())
    sim = test_convergence(dts, prob, Alg, reltol = 1.0e-16)
    @test sim.𝒪est[:l2] ≈ alg_order(Alg) atol = 0.5

    dts = 1 ./ 2 .^ (8:-1:4)
    sim = test_convergence(dts, prob, ETDRK4(), dense_errors = true)
    @test sim.𝒪est[:l2] ≈ 4 atol = 0.1
    @test sim.𝒪est[:L2] ≈ 4 atol = 0.1
end

@info "CFNLIRK3() is broken"

@testset "EPIRK Out-of-place" begin
    println("EPIRK Out-of-place")
    # Setup nonlinear problem
    A = [-2.0 1.0; 1.0 -2.0]
    f = (u, p, t) -> A * u - u .^ 3
    jac = (u, p, t) -> A - [3u[1]^2 0.0; 0.0 3u[2]^2]
    fun = ODEFunction(f; jac)
    Random.seed!(0)
    u0 = rand(2)
    tspan = (0.0, 1.0)
    prob = ODEProblem(fun, u0, tspan)
    # Setup approximate solution
    test_setup = Dict(:alg => Vern9(), :reltol => 1.0e-16, :abstol => 1.0e-16)
    # Convergence simulation
    dts = 1 ./ 2 .^ (7:-1:4)
    Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2]
    for Alg in Algs
        sim = analyticless_test_convergence(
            dts, prob, Alg(adaptive_krylov = false),
            test_setup
        )
        if Alg == EPIRK5s3
            @test_broken sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
        else
            @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
        end
    end
end

@testset "EPIRK Inplace" begin
    println("EPIRK Inplace")
    # Setup nonlinear problem
    A = [-2.0 1.0; 1.0 -2.0]
    f = (du, u, p, t) -> (mul!(du, A, u); du .-= u .^ 3)
    jac_update! = (J, u, p, t) -> (copyto!(J, A); J[1, 1] -= 3u[1]^2; J[2, 2] -= 3u[2]^2)
    jac_prototype = MatrixOperator(zeros(2, 2); update_func! = jac_update!)
    fun = ODEFunction(f; jac_prototype)
    Random.seed!(0)
    u0 = rand(2)
    tspan = (0.0, 1.0)
    prob = ODEProblem(fun, u0, tspan)
    # Setup approximate solution
    test_setup = Dict(:alg => Vern9(), :reltol => 1.0e-16, :abstol => 1.0e-16)
    # Convergence simulation
    dts = 1 ./ 2 .^ (7:-1:4)
    Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2]
    for Alg in Algs
        sim = analyticless_test_convergence(
            dts, prob, Alg(adaptive_krylov = false),
            test_setup
        )
        if Alg == EPIRK5s3
            @test_broken sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
        else
            @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
        end
    end
end

@testset "ExpRK Autodiff" begin
    println("ExpRK Autodiff")
    # Setup nonlinear problem (without explicit jacobian)
    A = [-2.0 1.0; 1.0 -2.0]
    f = (u, p, t) -> A * u - u .^ 3
    f_ip = (du, u, p, t) -> (mul!(du, A, u); du .-= u .^ 3)
    fun = ODEFunction(f)
    fun_ip = ODEFunction(f_ip)
    Random.seed!(0)
    u0 = rand(2)
    tspan = (0.0, 1.0)
    prob = ODEProblem(fun, u0, tspan)
    prob_ip = ODEProblem(fun_ip, u0, tspan)
    # Setup approximate solution
    test_setup = Dict(:alg => Vern9(), :reltol => 1.0e-16, :abstol => 1.0e-16)
    # Convergence simulation
    dts = 1 ./ 2 .^ (7:-1:4)
    sim = analyticless_test_convergence(dts, prob, HochOst4(krylov = true), test_setup)
    @test sim.𝒪est[:l2] ≈ 4 atol = 0.1
    sim = analyticless_test_convergence(dts, prob_ip, HochOst4(krylov = true), test_setup)
    @test sim.𝒪est[:l2] ≈ 4 atol = 0.1
    sim = analyticless_test_convergence(
        dts, prob, EPIRK5P1(adaptive_krylov = false),
        test_setup
    )
    @test sim.𝒪est[:l2] ≈ 5 atol = 0.1
    sim = analyticless_test_convergence(
        dts, prob_ip, EPIRK5P1(adaptive_krylov = false),
        test_setup
    )
    @test sim.𝒪est[:l2] ≈ 5 atol = 0.1
end

@testset "Adaptive Exprb Out-of-place" begin
    println("Adaptive Exprb Out-of-place")
    # Setup nonlinear problem
    A = [-2.0 1.0; 1.0 -2.0]
    f = (u, p, t) -> A * u - u .^ 3
    jac = (u, p, t) -> A - [3u[1]^2 0.0; 0.0 3u[2]^2]
    fun = ODEFunction(f; jac)
    Random.seed!(0)
    u0 = rand(2)
    tspan = (0.0, 1.0)
    prob = ODEProblem(fun, u0, tspan)
    # Setup approximate solution
    test_setup = Dict(:alg => Vern9(), :reltol => 1.0e-16, :abstol => 1.0e-16)
    # Convergence simulation
    dts = 1 ./ 2 .^ (7:-1:4)
    Algs = [Exprb32, Exprb43]
    for Alg in Algs
        sim = analyticless_test_convergence(dts, prob, Alg(), test_setup)
        @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
    end
end

@testset "Adaptive Exprb Inplace" begin
    println("Adaptive Exprb Inplace")
    # Setup nonlinear problem
    A = [-2.0 1.0; 1.0 -2.0]
    f = (du, u, p, t) -> (mul!(du, A, u); du .-= u .^ 3)
    jac_update! = (J, u, p, t) -> (copyto!(J, A); J[1, 1] -= 3u[1]^2; J[2, 2] -= 3u[2]^2)
    jac_prototype = MatrixOperator(zeros(2, 2); update_func! = jac_update!)
    fun = ODEFunction(f; jac_prototype)
    Random.seed!(0)
    u0 = rand(2)
    tspan = (0.0, 1.0)
    prob = ODEProblem(fun, u0, tspan)
    # Setup approximate solution
    test_setup = Dict(:alg => Vern9(), :reltol => 1.0e-16, :abstol => 1.0e-16)
    # Convergence simulation
    dts = 1 ./ 2 .^ (7:-1:4)
    Algs = [Exprb32, Exprb43]
    for Alg in Algs
        sim = analyticless_test_convergence(dts, prob, Alg(), test_setup)
        @test sim.𝒪est[:l2] ≈ alg_order(Alg()) atol = 0.1
    end
end

@testset "Exprb methods with the default in-place problem specialization" begin
    f! = (du, u, p, t) -> (du .= -u .+ u .^ 2 ./ 10)
    prob = ODEProblem(f!, [1.0, 0.5], (0.0, 1.0))
    prob_full = ODEProblem{true, SciMLBase.FullSpecialize}(f!, [1.0, 0.5], (0.0, 1.0))
    for alg in (Exprb32(), Exprb43())
        sol = solve(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        sol_full = solve(prob_full, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end] ≈ sol_full.u[end]
    end
end

@testset "Embedded ExpRK pair tableaus" begin
    # Both estimates of a pair have row sum phi_1(hA), so the weights of their
    # difference sum to zero, and each internal stage row sums to
    # c_i * phi_1(c_i * hA).
    A = [-3.0 1.0; 1.0 -4.0]
    dt = 0.3
    p1(c) = phi(c * dt * A, 1)[2]
    Z = zeros(2, 2)

    A21, A31, A32, B1, B2, B3, E1, E2, E3, E4 = expRK_operators(ERK32ZB(), dt, A)
    @test A21 ≈ (1 / 2) * p1(1 / 2)
    @test A31 + A32 ≈ (3 / 4) * p1(3 / 4)
    @test B1 + B2 + B3 ≈ p1(1)
    @test E1 + E2 + E3 + E4 ≈ Z atol = 1.0e-12

    ops = expRK_operators(ERK43ZB(), dt, A)
    A21, A31, A32, A41, A42, A43, A51, A52, A53, A54 = ops[1:10]
    B1, B2, B3, B4, B5 = ops[11:15]
    E1, E2, E3, E4, E5 = ops[16:20]
    @test A21 ≈ (1 / 6) * p1(1 / 6)
    @test A31 + A32 ≈ (1 / 2) * p1(1 / 2)
    @test A41 + A42 + A43 ≈ (1 / 2) * p1(1 / 2)
    @test A51 + A52 + A53 + A54 ≈ p1(1)
    @test B1 + B2 + B3 + B4 + B5 ≈ p1(1)
    @test E1 + E2 + E3 + E4 + E5 ≈ Z atol = 1.0e-12
end

@testset "Embedded ExpRK pairs adaptive" begin
    println("Embedded ExpRK pairs adaptive")
    μ = 1.01
    A = [2.0 -1.0; -1.0 2.0]
    u0 = [0.3, -0.2]
    analytic = (u0, p, t) -> exp((A + μ * I) * t) * u0
    prob = SplitODEProblem(
        SplitFunction(
            MatrixOperator(A), (du, u, p, t) -> du .= μ .* u; analytic
        ), u0, (0.0, 1.0)
    )
    exact = analytic(u0, nothing, 1.0)

    for Alg in (ERK32ZB, ERK43ZB)
        loose = solve(prob, Alg(); abstol = 1.0e-4, reltol = 1.0e-4)
        tight = solve(prob, Alg(); abstol = 1.0e-9, reltol = 1.0e-9)
        @test SciMLBase.successful_retcode(tight)
        @test norm(tight.u[end] - exact) < norm(loose.u[end] - exact)
        @test norm(tight.u[end] - exact) < 1.0e-6
        # the operators are rebuilt whenever dt moves, so a step size that really
        # varies is what exercises that path
        interior = diff(tight.t)[1:(end - 1)]
        @test maximum(interior) / minimum(interior) > 1.5
        # Krylov and cached operators are two routes to the same step
        kry = solve(prob, Alg(krylov = true); abstol = 1.0e-9, reltol = 1.0e-9)
        @test kry.t ≈ tight.t
        @test kry.u[end] ≈ tight.u[end]
    end
end

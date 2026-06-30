# This definitely needs cleaning
using OrdinaryDiffEqBDF, ODEProblemLibrary, DiffEqDevTools, ADTypes, LinearSolve, StaticArrays
using OrdinaryDiffEqNonlinearSolve: NLFunctional, NLAnderson, NonlinearSolveAlg
using Test, Random
Random.seed!(100)

testTol = 0.2
dts = 1 .// 2 .^ (9:-1:5)
dts3 = 1 .// 2 .^ (12:-1:7)

if isempty(VERSION.prerelease)
    using Enzyme
end

@testset "Implicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in
    1:2

    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    sim = test_convergence(dts, prob, ABDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    sim = test_convergence(dts, prob, ABDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    # QBDF
    sim = test_convergence(dts, prob, QBDF1())
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    sim = test_convergence(dts, prob, QBDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    # QNDF
    sim = test_convergence(dts, prob, QNDF1())
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    sim = test_convergence(dts, prob, QNDF1(autodiff = AutoFiniteDiff()))
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    @test sim.𝒪est[:l2] ≈ 1 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            QNDF1(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                )
            )
        )
        @test sim.𝒪est[:final] ≈ 1 atol = testTol
        @test sim.𝒪est[:l2] ≈ 1 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol

        # Note: tight reltol/abstol (1e-14) is unusable for BDF-1 because the
        # method's local error floor exceeds machine precision, causing Newton
        # ConvergenceFailure on all but the finest dt. Leave reltol/abstol at the
        # defaults; without the LinearSolve residual-tolerance coupling described in
        # #3373, QNDF1's observed order already lands near 1 on prob_ode_2Dlinear.
        sim = test_convergence(
            dts,
            prob,
            QNDF1(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 1 atol = testTol
        @test sim.𝒪est[:l2] ≈ 1 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 1 atol = testTol
    end

    sim = test_convergence(dts3, prob, QNDF2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol

    sim = test_convergence(dts, prob, QNDF2(nlsolve = NLFunctional()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:l2] ≈ 2 atol = testTol
    @test sim.𝒪est[:l∞] ≈ 2 atol = testTol
    @test_nowarn solve(prob, QNDF())

    # MEBDF2
    sim21 = test_convergence(dts, prob, MEBDF2(extrapolant = :linear))
    @test sim21.𝒪est[:final] ≈ 2 atol = testTol

    sim22 = test_convergence(dts, prob, MEBDF2(nlsolve = NLFunctional()), reltol = 1.0e-2)
    @test sim22.𝒪est[:final] ≈ 2 atol = testTol

    sim23 = test_convergence(dts, prob, MEBDF2(nlsolve = NLAnderson()), reltol = 1.0e-2)
    @test sim23.𝒪est[:final] ≈ 2 atol = testTol

    sim24 = test_convergence(
        dts, prob, MEBDF2(nlsolve = NonlinearSolveAlg()), reltol = 1.0e-2
    )
    @test sim24.𝒪est[:final] ≈ 2 atol = testTol

    #FBDF
    @test_nowarn solve(prob, FBDF())

    # MOOSE234 — VSVO method, test_convergence with fixed dt is not applicable
    @test_nowarn solve(prob, MOOSE234())

    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
end

@testset "MOOSE234 accuracy (scalar exponential)" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    exact = exp(-1.0)

    sol_mod = solve(prob, MOOSE234(), abstol = 1e-6, reltol = 1e-6)
    @test isapprox(sol_mod[end], exact, rtol = 1e-4)

    sol_tight = solve(prob, MOOSE234(), abstol = 1e-10, reltol = 1e-10)
    @test isapprox(sol_tight[end], exact, rtol = 1e-6)

    @test abs(sol_tight[end] - exact) < abs(sol_mod[end] - exact)
end

@testset "MOOSE234 accuracy (in-place 2D)" begin
    fiip = (du, u, p, t) -> du .= p .* u
    u0 = [1.0, 2.0]
    prob = ODEProblem(fiip, u0, (0.0, 1.0), -1.0)

    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
    @test isapprox(sol.u[end], exp(-1.0) .* u0, rtol = 1e-2)
end

@testset "MOOSE234 vs FBDF accuracy" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    exact = exp(-1.0)

    sol_moose = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    sol_fbdf = solve(prob, FBDF(), abstol = 1e-8, reltol = 1e-8)

    err_moose = abs(sol_moose[end] - exact)
    err_fbdf = abs(sol_fbdf[end] - exact)

    @test err_moose < 1e-6
    @test err_fbdf < 1e-6
end

@testset "Static Array (SVector) Tests" begin
    f_oop(u, p, t) = -0.5 * u
    u0_sv = SVector(1.0, 2.0)
    prob_sv = ODEProblem(f_oop, u0_sv, (0.0, 1.0))

    for alg in (QNDF(), QNDF1(), QNDF2(), FBDF(), MOOSE234())
        name = nameof(typeof(alg))
        @testset "$name" begin
            sol = solve(prob_sv, alg, abstol = 1.0e-8, reltol = 1.0e-8)
            @test sol.u[end] isa SVector
            @test isapprox(sol.u[end], exp(-0.5) * u0_sv, rtol = 1.0e-3)
        end
    end

    # Scalar
    prob_scalar = ODEProblem(f_oop, 1.0, (0.0, 1.0))
    for alg in (QNDF(), FBDF(), MOOSE234())
        name = nameof(typeof(alg))
        @testset "$name scalar" begin
            sol = solve(prob_scalar, alg, abstol = 1.0e-8, reltol = 1.0e-8)
            @test sol.u[end] isa Number
            @test isapprox(sol.u[end], exp(-0.5), rtol = 1.0e-5)
        end
    end
end

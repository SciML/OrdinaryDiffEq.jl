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
end

@testset "SBDF" begin
    # Regression for #4329: SBDF2/3/4 silently ran below their nominal order.
    # `measured_order` runs a fixed-dt, non-adaptive solve at two step sizes
    # and infers the order from the error ratio (test_convergence's
    # simulation machinery assumes a plain ODEProblem; SBDF is IMEX-only).
    # Errors are reduced with `maximum(abs, ...)` so this works for both the
    # scalar and vector (in-place) problems below.
    exact(t) = exp(-t)
    function measured_order(prob, alg, dts)
        errs = [
            maximum(abs, solve(prob, alg; dt, adaptive = false).u[end] .- exact(prob.tspan[2]))
                for dt in dts
        ]
        return log(errs[1] / errs[2]) / log(dts[1] / dts[2])
    end

    @testset "pure-BDF equivalent (f2 ≡ 0), out-of-place" begin
        prob = SplitODEProblem((u, p, t) -> -u, (u, p, t) -> 0.0, 1.0, (0.0, 1.0))
        for (alg, order) in ((SBDF2(), 2), (SBDF3(), 2), (SBDF4(), 2))
            # Capped at 2 regardless of alg.order: the first alg.order - 1
            # steps bootstrap through the lower-order family members (no
            # trusted multi-point history exists yet), and the first step's
            # O(dt^2) local error survives in the global error for every
            # later step -- this is expected multistep start-up behavior,
            # not a defect (see #4329's "related observations").
            @test measured_order(prob, alg, (1.0e-2, 1.0e-3)) ≈ order atol = 0.2
        end
    end

    @testset "genuine IMEX split, out-of-place and in-place" begin
        # The first step happens to be the trapezoidal rule on an equal
        # split, so the start-up cap here is 3, not 2.
        prob = SplitODEProblem((u, p, t) -> -0.5 * u, (u, p, t) -> -0.5 * u, 1.0, (0.0, 1.0))
        prob_ip = SplitODEProblem(
            (du, u, p, t) -> (du .= -0.5 .* u; nothing),
            (du, u, p, t) -> (du .= -0.5 .* u; nothing),
            [1.0], (0.0, 1.0)
        )
        for (p, order) in ((prob, 2), (prob_ip, 2))
            @test measured_order(p, SBDF2(), (1.0e-2, 1.0e-3)) ≈ order atol = 0.2
        end
        for p in (prob, prob_ip)
            @test measured_order(p, SBDF3(), (1.0e-2, 1.0e-3)) ≈ 3 atol = 0.2
            @test measured_order(p, SBDF4(), (1.0e-2, 1.0e-3)) ≈ 3 atol = 0.2
        end
    end

    @testset "step-size change does not corrupt the fixed-coefficient history" begin
        # A clipped final step (dt landing just short of tspan[2]) must not
        # pollute the result -- this reproduces the ~4200x error inflation
        # from #4329 (defect 3) if unfixed.
        prob = SplitODEProblem((u, p, t) -> -u, (u, p, t) -> 0.0, 1.0, (0.0, 2.0))
        err_clipped = abs(
            solve(prob, SBDF2(); dt = 1.0e-3, adaptive = false).u[end] - exact(2.0)
        )
        err_exact_binary = abs(
            solve(prob, SBDF2(); dt = 2.0^-10, adaptive = false).u[end] - exact(2.0)
        )
        @test err_clipped < 5 * err_exact_binary
    end
end

@testset "Static Array (SVector) Tests" begin
    f_oop(u, p, t) = -0.5 * u
    u0_sv = SVector(1.0, 2.0)
    prob_sv = ODEProblem(f_oop, u0_sv, (0.0, 1.0))

    for alg in (QNDF(), QNDF1(), QNDF2(), FBDF())
        name = nameof(typeof(alg))
        @testset "$name" begin
            sol = solve(prob_sv, alg, abstol = 1.0e-8, reltol = 1.0e-8)
            @test sol.u[end] isa SVector
            @test isapprox(sol.u[end], exp(-0.5) * u0_sv, rtol = 1.0e-3)
        end
    end

    # Scalar
    prob_scalar = ODEProblem(f_oop, 1.0, (0.0, 1.0))
    for alg in (QNDF(), FBDF())
        name = nameof(typeof(alg))
        @testset "$name scalar" begin
            sol = solve(prob_scalar, alg, abstol = 1.0e-8, reltol = 1.0e-8)
            @test sol.u[end] isa Number
            @test isapprox(sol.u[end], exp(-0.5), rtol = 1.0e-5)
        end
    end
end

using GlobalDiffEq, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, LinearAlgebra
using Test
import SciMLBase

lv!(du, u, p, t) = (du[1] = 1.5u[1] - u[1] * u[2]; du[2] = -3.0u[2] + u[1] * u[2]; nothing)
lv(u, p, t) = [1.5u[1] - u[1] * u[2], -3.0u[2] + u[1] * u[2]]
const lv_tspan = (0.0, 10.0)

const EQUATIONS = (DefectCorrection, ErrorTransport)
const MODES = (InterpolatingMode, SimultaneousMode)

@testset "Companion global error estimators" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    prob_oop = ODEProblem(lv, [1.0, 1.0], lv_tspan)
    ref = solve(prob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)

    for equation in EQUATIONS, mode in MODES, p in (prob, prob_oop)
        for tol in (1.0e-4, 1.0e-6)
            alg = GlobalErrorEstimation(Tsit5(); equation, mode)
            est = global_error_estimate(p, alg; abstol = tol, reltol = tol)
            sol = solve(p, Tsit5(); abstol = tol, reltol = tol)
            true_err = norm(sol.u[end] - ref.u[end])
            @test est / true_err ≈ 1 rtol = 0.1
        end
    end
end

@testset "SimultaneousMode matches InterpolatingMode" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    for equation in EQUATIONS
        for tol in (1.0e-4, 1.0e-6)
            interp = global_error_estimate(
                prob, GlobalErrorEstimation(Tsit5(); equation, mode = InterpolatingMode);
                abstol = tol, reltol = tol
            )
            stream = global_error_estimate(
                prob, GlobalErrorEstimation(Tsit5(); equation, mode = SimultaneousMode);
                abstol = tol, reltol = tol
            )
            # Same companion ODE; the interpolating mode integrates it in one solve
            # over the whole span while the streaming mode restarts the companion
            # solver at every forward step, so they agree to a few 1e-3 (negligible
            # next to the ~10% accuracy of the estimate itself), not to machine tol.
            @test stream ≈ interp rtol = 5.0e-3
        end
    end
end

@testset "Companion global error control" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    ref = solve(prob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    gtol = 1.0e-7

    for equation in EQUATIONS, mode in MODES
        alg = GlobalErrorEstimation(Tsit5(); equation, mode, gtol)
        sol = solve(prob, alg; abstol = 1.0e-3, reltol = 1.0e-3)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] - ref.u[end]) <= gtol
    end
end

@testset "Error trajectory: $mode" for mode in MODES
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    ref = solve(prob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    alg = GlobalErrorEstimation(Tsit5(); mode)
    @test SciMLBase.has_global_error(alg)
    @test !SciMLBase.has_global_error(GlobalErrorEstimation(Tsit5(); equation = ErrorTransport))

    # Without gtol, the solve gives the error of each saved state at the given tolerances.
    sol = solve(prob, alg; abstol = 1.0e-6, reltol = 1.0e-6)
    @test length(sol.global_error) == length(sol.u)
    @test global_error_estimate(sol) === sol.global_error
    @test norm(sol.global_error[end]) ≈
        global_error_estimate(prob, alg; abstol = 1.0e-6, reltol = 1.0e-6)
    actual = [ref(t) .- u for (t, u) in zip(sol.t, sol.u)]
    @test all(
        norm(e - a) <= 1.0e-3 * maximum(norm, actual) for (e, a) in zip(sol.global_error, actual)
    )

    # Under global error control, the final solution has the error trajectory too.
    gsol = solve(
        prob, GlobalErrorEstimation(Tsit5(); mode, gtol = 1.0e-7); abstol = 1.0e-3, reltol = 1.0e-3
    )
    @test length(gsol.global_error) == length(gsol.u)
    gactual = [ref(t) .- u for (t, u) in zip(gsol.t, gsol.u)]
    @test all(
        norm(e - a) <= 1.0e-3 * maximum(norm, gactual) for (e, a) in zip(gsol.global_error, gactual)
    )
end

@testset "Stiff problem with an implicit solver" begin
    # Van der Pol with mu = 1000. At the fast transition near t = 807 the global error
    # passes through a peak 1e7 times its final value, so a loose default companion
    # reltol gives a wrong estimate.
    svdp!(du, u, p, t) = (du[1] = u[2]; du[2] = 1.0e3 * ((1 - u[1]^2) * u[2]) - u[1]; nothing)
    prob = ODEProblem(svdp!, [2.0, 0.0], (0.0, 1000.0))
    ref = solve(prob, Rodas5P(); abstol = 1.0e-12, reltol = 1.0e-12)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-4, reltol = 1.0e-4)
    gtol = 1.0e-6
    for mode in MODES
        est = global_error_estimate(
            prob, GlobalErrorEstimation(Rodas5P(); mode); abstol = 1.0e-4, reltol = 1.0e-4
        )
        @test est / norm(sol.u[end] - ref.u[end]) ≈ 1 rtol = 1.0e-2
        # The outer solve wraps f for the ForwardDiff traits of the inner solver, which
        # differentiates f for its Jacobian.
        gsol = solve(
            prob, GlobalErrorEstimation(Rodas5P(); mode, gtol); abstol = 1.0e-3, reltol = 1.0e-3
        )
        @test SciMLBase.successful_retcode(gsol)
        @test norm(gsol.u[end] - ref.u[end]) <= gtol
    end
end

@testset "Companion estimator argument validation" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    alg = GlobalErrorEstimation(Tsit5())

    @test !SciMLBase.allows_arbitrary_number_types(alg)
    @test !SciMLBase.allowscomplex(alg)
    @test !SciMLBase.isautodifferentiable(alg)
    @test alg.equation === DefectCorrection
    @test alg.mode === InterpolatingMode
    @test GlobalErrorEstimation(Tsit5(); equation = ErrorTransport).equation === ErrorTransport
    @test GlobalErrorEstimation(Tsit5(); mode = SimultaneousMode).mode === SimultaneousMode
    @test_throws ArgumentError GlobalErrorEstimation(Tsit5(); gtol = -1.0)
    @test_throws ArgumentError GlobalErrorEstimation(Tsit5(); maxiters = 0)
    @test_throws ArgumentError GlobalErrorEstimation(Tsit5(); safety = 2.0)
    # Without gtol, ErrorTransport has nothing to return (DefectCorrection returns the
    # error trajectory).
    @test_throws ArgumentError solve(prob, GlobalErrorEstimation(Tsit5(); equation = ErrorTransport))

    # Callbacks are supported only by DefectCorrection (see companion_events_tests.jl).
    callback = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    callback_prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan; callback)
    for mode in MODES
        @test global_error_estimate(callback_prob, GlobalErrorEstimation(Tsit5(); mode)) > 0
        unsupported = GlobalErrorEstimation(Tsit5(); equation = ErrorTransport, mode)
        @test_throws ArgumentError global_error_estimate(callback_prob, unsupported)
        @test_throws ArgumentError global_error_estimate(prob, unsupported; callback)
    end

    backwards_prob = ODEProblem(lv!, [1.0, 1.0], (10.0, 0.0))
    @test_throws ArgumentError global_error_estimate(backwards_prob, alg)

    # Both integration modes require a solver with genuine dense output; a solution
    # left with the trivial linear-interpolation fallback is rejected rather than
    # silently differentiated into a wrong defect.
    @testset "dense-output requirement" begin
        t = [0.0, 0.5, 1.0]
        u = [[1.0], [0.61], [0.37]]
        nondense = SciMLBase.build_solution(
            prob, Tsit5(), t, u;
            interp = SciMLBase.LinearInterpolation(t, u), dense = false,
            retcode = SciMLBase.ReturnCode.Success
        )
        @test_throws ArgumentError GlobalDiffEq._validate_estimation_solution(
            nondense, "GlobalErrorEstimation"
        )
        @test_throws ArgumentError GlobalDiffEq._validate_streaming_interpolation(
            nondense, "GlobalErrorEstimation"
        )
    end
end

@testset "Companion solver warning" begin
    # The companion runs at tight tolerances, which a low-order default solver meets poorly.
    @test_logs (:warn, r"solves the companion with Rosenbrock32") GlobalErrorEstimation(Rosenbrock32())
    @test_logs GlobalErrorEstimation(Rosenbrock32(); companion_alg = Tsit5())
    @test_logs GlobalErrorEstimation(Tsit5())
    @test GlobalErrorEstimation(Tsit5()).companion_alg isa Tsit5
end

using GlobalDiffEq, OrdinaryDiffEqTsit5, LinearAlgebra
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
    # gtol is required to solve with the wrapper
    @test_throws ArgumentError solve(prob, alg)

    callback = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    callback_prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan; callback)
    @test_throws ArgumentError global_error_estimate(callback_prob, alg)

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

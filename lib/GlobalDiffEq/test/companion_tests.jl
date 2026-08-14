using GlobalDiffEq, OrdinaryDiffEqTsit5, LinearAlgebra
using Test
import SciMLBase

lv!(du, u, p, t) = (du[1] = 1.5u[1] - u[1] * u[2]; du[2] = -3.0u[2] + u[1] * u[2]; nothing)
lv(u, p, t) = [1.5u[1] - u[1] * u[2], -3.0u[2] + u[1] * u[2]]
const lv_tspan = (0.0, 10.0)

@testset "Companion global error estimators" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    ref = solve(prob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)

    for alg in (GlobalErrorTransport(Tsit5()),)
        for tol in (1.0e-4, 1.0e-6)
            est = global_error_estimate(prob, alg; abstol = tol, reltol = tol)
            sol = solve(prob, Tsit5(); abstol = tol, reltol = tol)
            true_err = norm(sol.u[end] - ref.u[end])
            @test est / true_err ≈ 1 rtol = 0.1
        end
    end

    # out-of-place problems
    prob_oop = ODEProblem(lv, [1.0, 1.0], lv_tspan)
    for alg in (GlobalErrorTransport(Tsit5()),)
        est = global_error_estimate(prob_oop, alg; abstol = 1.0e-6, reltol = 1.0e-6)
        sol = solve(prob_oop, Tsit5(); abstol = 1.0e-6, reltol = 1.0e-6)
        true_err = norm(sol.u[end] - ref.u[end])
        @test est / true_err ≈ 1 rtol = 0.1
    end
end

@testset "Companion global error control" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)
    ref = solve(prob, Tsit5(); abstol = 1.0e-13, reltol = 1.0e-13)
    gtol = 1.0e-7

    for alg in (
            GlobalErrorTransport(Tsit5(); gtol),
        )
        sol = solve(prob, alg; abstol = 1.0e-3, reltol = 1.0e-3)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] - ref.u[end]) <= gtol
    end
end

@testset "Companion estimator argument validation" begin
    prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan)

    for Alg in (GlobalErrorTransport,)
        alg = Alg(Tsit5())
        @test !SciMLBase.allows_arbitrary_number_types(alg)
        @test !SciMLBase.allowscomplex(alg)
        @test !SciMLBase.isautodifferentiable(alg)
        @test_throws ArgumentError Alg(Tsit5(); gtol = -1.0)
        @test_throws ArgumentError Alg(Tsit5(); maxiters = 0)
        @test_throws ArgumentError Alg(Tsit5(); safety = 2.0)
        # gtol is required to solve with the wrapper
        @test_throws ArgumentError solve(prob, alg)

        callback = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
        callback_prob = ODEProblem(lv!, [1.0, 1.0], lv_tspan; callback)
        @test_throws ArgumentError global_error_estimate(callback_prob, alg)

        backwards_prob = ODEProblem(lv!, [1.0, 1.0], (10.0, 0.0))
        @test_throws ArgumentError global_error_estimate(backwards_prob, alg)
    end
end

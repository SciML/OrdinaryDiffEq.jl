using GlobalDiffEq, OrdinaryDiffEqTsit5, Random
using SciMLSensitivity, QuadGK
using Test
import SciMLBase

@testset "Adjoint error estimation and control" begin
    function linear!(du, u, p, t)
        du[1] = p * u[1]
        return nothing
    end

    rate = 2.0
    tspan = (0.0, 2.0)
    linear_prob = ODEProblem(linear!, [1.0], tspan, rate)
    estimator_alg = GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(123))

    coarse_sol = solve(
        linear_prob, Tsit5();
        abstol = 1.0e-4, reltol = 1.0e-4, dense = true, save_everystep = true
    )
    estimated_error = adjoint_error_estimate(
        linear_prob, estimator_alg; abstol = 1.0e-4, reltol = 1.0e-4
    )
    actual_error = abs(coarse_sol.u[end][1] - exp(rate * tspan[2]))

    @test estimated_error > 0
    @test estimated_error / actual_error ≈ 1 rtol = 0.15

    gtol = 1.0e-6
    controlled_sol = solve(
        linear_prob,
        GlobalAdjoint(Tsit5(); gtol, samples = 1, rng = Xoshiro(321));
        abstol = 1.0e-3, reltol = 1.0e-3
    )
    @test abs(controlled_sol.u[end][1] - exp(rate * tspan[2])) <= gtol
    @test controlled_sol.prob.p === rate

    forcing = t -> 0.1sin(t)
    forced(u, p, t) = [-u[1] + p(t)]

    forced_prob = ODEProblem(forced, [0.0], (0.0, 1.0), forcing)
    callable_parameter_estimate = adjoint_error_estimate(
        forced_prob, GlobalAdjoint(Tsit5(); samples = 1, rng = Xoshiro(456));
        abstol = 1.0e-4, reltol = 1.0e-4
    )
    @test isfinite(callable_parameter_estimate)
    @test callable_parameter_estimate >= 0

    adjoint_alg = GlobalAdjoint(Tsit5())
    @test !SciMLBase.allows_arbitrary_number_types(adjoint_alg)
    @test !SciMLBase.allowscomplex(adjoint_alg)
    @test !SciMLBase.isautodifferentiable(adjoint_alg)
    @test_throws ArgumentError GlobalAdjoint(Tsit5(); samples = 0)
    @test_throws ArgumentError GlobalAdjoint(Tsit5(); gtol = -1.0)
    @test_throws ArgumentError GlobalAdjoint(Tsit5(); sensealg = :not_a_sensealg)
    @test_throws ArgumentError solve(linear_prob, adjoint_alg)

    callback = DiscreteCallback((u, t, integrator) -> false, integrator -> nothing)
    callback_prob = ODEProblem(linear!, [1.0], tspan, rate; callback)
    @test_throws ArgumentError adjoint_error_estimate(callback_prob, adjoint_alg)
end

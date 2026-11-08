using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, LinearSolve, Test
using LinearAlgebra: lu!
using SciMLBase: SciMLBase

# The Jacobian is NaN everywhere, so the Newton matrix W is non-finite on every attempt. LAPACK's LU throws on non-finite input; the Newton solve
# must report divergence instead, and the solve must stop with a failure retcode rather
# than an exception.
f!(du, u, p, t) = (du .= -u; nothing)
jac!(J, u, p, t) = (fill!(J, NaN); nothing)
prob = ODEProblem(ODEFunction(f!; jac = jac!), [1.0, 1.0], (0.0, 1.0))
algs = (ImplicitEuler, TRBDF2, KenCarp4, FBDF, QNDF)
linsolves = (LUFactorization(), GenericFactorization(lu!))

@testset "Non-finite W fails the Newton solve" begin
    @testset "$(nameof(alg)) with $(nameof(typeof(linsolve)))" for alg in algs,
            linsolve in linsolves

        sol = solve(prob, alg(; linsolve))
        @test !SciMLBase.successful_retcode(sol)
        @test sol.t[end] == 0
        @test sol.u[end] == prob.u0
    end
end

@testset "Finite W with generic factorization" begin
    finite_prob = ODEProblem(f!, copy(prob.u0), prob.tspan)
    for alg in algs
        reference = solve(
            finite_prob, alg(; linsolve = LUFactorization()); adaptive = false, dt = 0.1
        )
        sol = solve(
            finite_prob, alg(; linsolve = GenericFactorization(lu!)); adaptive = false,
            dt = 0.1
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == reference.t
        @test all(
            isapprox(u, v; rtol = 1.0e-12, atol = 1.0e-14)
                for (u, v) in zip(sol.u, reference.u)
        )
    end
end

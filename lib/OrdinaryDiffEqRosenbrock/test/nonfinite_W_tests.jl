using OrdinaryDiffEqRosenbrock, LinearSolve, Test
using LinearAlgebra: lu!
using SciMLBase: SciMLBase

# The Jacobian is NaN once the decaying solution drops below 1/2, so every step from such
# a point forms a non-finite W. LAPACK's LU throws on non-finite input; the step must be
# rejected instead, and the solve must stop with a failure retcode rather than an exception.
f!(du, u, p, t) = (du .= -u; nothing)
jac!(J, u, p, t) = (J .= u[1] < 0.5 ? NaN : -1.0; nothing)
prob = ODEProblem(ODEFunction(f!; jac = jac!), [1.0, 1.0], (0.0, 2.0))
algs = (Rodas5P, Rosenbrock23, Rosenbrock32)
linsolves = (LUFactorization(), GenericFactorization(lu!))

@testset "Non-finite W rejects the step" begin
    @testset "$(nameof(alg)) with $(nameof(typeof(linsolve)))" for alg in algs,
            linsolve in linsolves

        sol = solve(prob, alg(; linsolve))
        @test !SciMLBase.successful_retcode(sol)
        @test log(2) < sol.t[end] < 1.0
        @test all(isfinite, sol.u[end])
    end
end

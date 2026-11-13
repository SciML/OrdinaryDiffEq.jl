using OrdinaryDiffEqRosenbrock, LinearSolve, Test
using LinearAlgebra: lu!, I
using SciMLBase: SciMLBase

# The Jacobian is NaN once the decaying solution drops below 1/2, so every step from such
# a point forms a non-finite W. LAPACK's LU throws on non-finite input; the step must fail
# instead, and the solve must stop with a failure retcode rather than an exception. With
# fixed steps or `force_dtmin` an error-estimate rejection is not binding, so the failure
# must also stop those solves rather than advance time with an untouched state.
f!(du, u, p, t) = (du .= -u; nothing)
jac!(J, u, p, t) = (u[1] < 0.5 ? fill!(J, NaN) : copyto!(J, -I); nothing)
prob = ODEProblem(ODEFunction(f!; jac = jac!), [1.0, 1.0], (0.0, 2.0))
algs = (Rodas5P, Rosenbrock23, Rosenbrock32)
linsolves = (LUFactorization(), GenericFactorization(lu!))
modes = (
    adaptive = (;),
    fixed = (; adaptive = false, dt = 0.1),
    force_dtmin = (; dt = 0.1, dtmin = 0.1, force_dtmin = true),
)

@testset "Non-finite W fails the step" begin
    @testset "$(nameof(alg)) with $(nameof(typeof(linsolve))), $mode" for alg in algs,
            linsolve in linsolves, mode in keys(modes)

        sol = solve(prob, alg(; linsolve); modes[mode]...)
        @test !SciMLBase.successful_retcode(sol)
        @test log(2) < sol.t[end] < 1.0
        @test sol.u[end] ≈ fill(exp(-sol.t[end]), 2) rtol = 1.0e-2
    end
end

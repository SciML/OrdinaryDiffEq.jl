# Regression tests for SciML/OrdinaryDiffEq.jl#3933: implicit solvers with an
# AbstractSciMLOperator Jacobian (islinearfunction path) and a concrete-matrix
# linear solver silently returned garbage, because do_newJW never requested a W
# refactorization for linear operators (W = J - M/(γdt) changes with dt even
# when J is constant). Requires SciMLOperators >= 1.24.4 for the companion fix
# that rebuilds WOperator._concrete_form on convert for operator Jacobians.
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqDifferentiation, LinearSolve, SciMLOperators
using SciMLBase, LinearAlgebra, Test

@testset "stale W with linear operator f and concrete-A linsolve" begin
    A = [-2.0 1.0; 1.0 -2.0]
    B = [0.0 0.5; -0.5 0.0]
    u0 = [1.0, 0.5]
    tspan = (0.0, 1.0)
    uref = exp((A + B) * tspan[2]) * u0

    f2! = (du, u, p, t) -> mul!(du, B, u)
    prob_split = SplitODEProblem(MatrixOperator(A), f2!, u0, tspan)

    @testset "$(nameof(typeof(alg))) split, LU" for alg in (
            KenCarp3(linsolve = LUFactorization()),
            KenCarp4(linsolve = LUFactorization()),
            KenCarp5(linsolve = LUFactorization()),
        )
        sol = solve(prob_split, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] - uref) / norm(uref) < 1.0e-3
        # the broken path accepted the whole interval in ~2 giant steps
        @test length(sol.t) > 5
    end

    # Non-split linear-operator problem hits the same do_newJW islin branch
    prob_lin = ODEProblem(ODEFunction(MatrixOperator(A + B)), u0, tspan)
    sol = solve(prob_lin, TRBDF2(linsolve = LUFactorization()); abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test norm(sol.u[end] - uref) / norm(uref) < 1.0e-3

    # The Krylov path was always correct; make sure it stays that way
    sol_krylov = solve(
        prob_split, KenCarp4(linsolve = KrylovJL_GMRES());
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test norm(sol_krylov.u[end] - uref) / norm(uref) < 1.0e-3
end

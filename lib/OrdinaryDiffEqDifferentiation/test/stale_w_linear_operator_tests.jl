# Regression tests for SciML/OrdinaryDiffEq.jl#3933: implicit solvers with an
# AbstractSciMLOperator Jacobian (islinearfunction path) and a concrete-matrix
# linear solver silently returned garbage, because do_newJW never requested a W
# refactorization for linear operators (W = J - M/(γdt) changes with dt even
# when J is constant). Requires SciMLOperators >= 1.24.4 for the companion fix
# that rebuilds WOperator._concrete_form on convert for operator Jacobians.
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqDifferentiation, LinearSolve, SciMLOperators
using SciMLBase, LinearAlgebra, Test
using OrdinaryDiffEqSDIRK: NLNewton

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

@testset "fixed-step linear W is rebuilt whenever γdt changes (#4370)" begin
    A = [-20.0 1.0; 1.0 -20.0]
    u0 = [1.0, 0.5]
    prob = ODEProblem(ODEFunction(MatrixOperator(A)), u0, (0.0, 1.0))

    function stepped(cutoff)
        integ = init(
            prob, SDIRK2(nlsolve = NLNewton(new_W_dt_cutoff = cutoff));
            dt = 1.0e-2, adaptive = false
        )
        for k in 1:20
            set_proposed_dt!(integ, (1.0e-2, 1.15e-2)[mod1(k, 2)])
            step!(integ)
        end
        return integ.u, integ.t
    end

    u_exact_W, t = stepped(0.0)
    u_default, _ = stepped(0.2)
    @test u_default ≈ u_exact_W rtol = 1.0e-10

    uref = exp(A * t) * u0
    @test norm(u_default - uref) / norm(uref) < 0.05
end

@testset "split W (LHLFactorization) takes every γdt change (#4370)" begin
    A = [-20.0 1.0; 1.0 -20.0]
    u0 = [1.0, 0.5]
    prob = ODEProblem(ODEFunction(MatrixOperator(A)), u0, (0.0, 1.0))
    lhl(cutoff) = SDIRK2(linsolve = LHLFactorization(), nlsolve = NLNewton(new_W_dt_cutoff = cutoff))

    loose = solve(prob, lhl(0.2); abstol = 1.0e-6, reltol = 1.0e-6)
    exact = solve(prob, lhl(0.0); abstol = 1.0e-6, reltol = 1.0e-6)
    @test loose.u[end] == exact.u[end]
    @test loose.stats.naccept == exact.stats.naccept

    integ = init(prob, lhl(0.2); dt = 1.0e-2, adaptive = false)
    for k in 1:20
        set_proposed_dt!(integ, (1.0e-2, 1.15e-2)[mod1(k, 2)])
        step!(integ)
    end
    integ_exact = init(prob, lhl(0.0); dt = 1.0e-2, adaptive = false)
    for k in 1:20
        set_proposed_dt!(integ_exact, (1.0e-2, 1.15e-2)[mod1(k, 2)])
        step!(integ_exact)
    end
    @test integ.u == integ_exact.u
end

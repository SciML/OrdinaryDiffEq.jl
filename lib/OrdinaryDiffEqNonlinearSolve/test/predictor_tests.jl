using OrdinaryDiffEqNonlinearSolve:
    NLNewton, NLFunctional, NLAnderson, NonlinearSolveAlg, HomotopyNonlinearSolveAlg
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF
using NonlinearSolve: NewtonRaphson
using SciMLBase
using Test

# `nlsolve!` seeds each stage solve from the nlsolve algorithm's `predictor`
# callable: `predictor(uprev, p, t + c*dt, dt)` (out of place) guesses the stage
# *value*, and the solver maps it onto its own iterate `z`. These tests check
# that the seed actually reaches the stage residual's `ustep` for every stage of
# DIRK, IMEX, and coefficient-multistep methods, and for every nonlinear-solver
# backend, in place and out of place.

const SEED = 3.0

recording_f(evals) = (u, p, t) -> (u isa Float64 && push!(evals, u); u)

function recording_f!(evals)
    return (du, u, p, t) -> begin
        eltype(u) === Float64 && push!(evals, u[1])
        du .= u
        nothing
    end
end

@testset "callable stage predictor" begin
    @testset "seed reaches the residual (out of place)" begin
        for (make_alg, nsolves_per_step) in (
                (nl -> ImplicitEuler(nlsolve = nl), 1),
                (nl -> TRBDF2(nlsolve = nl), 2),
                (nl -> QNDF(nlsolve = nl), 1),
            )
            evals = Float64[]
            pred_ts = Float64[]
            pred = (uprev, p, t, dt) -> (push!(pred_ts, t); p)
            prob = ODEProblem(recording_f(evals), 1.0, (0.0, 1.0), SEED)
            sol = solve(
                prob, make_alg(NLNewton(predictor = pred));
                dt = 0.1, adaptive = false
            )
            @test sol.retcode == SciMLBase.ReturnCode.Success
            nsteps = length(sol.t) - 1
            @test nsteps == 10
            # the callable fires exactly once per stage solve, at the stage time
            @test length(pred_ts) == nsolves_per_step * nsteps
            @test issorted(pred_ts)
            @test all(0.0 .< pred_ts .<= 1.0)
            # and its output (it returned `p = SEED`) seeds the first residual eval
            @test count(==(SEED), evals) >= nsolves_per_step * nsteps
        end
    end

    @testset "seed reaches the residual (in place)" begin
        evals = Float64[]
        pred_ts = Float64[]
        function pred!(upred, uprev, p, t, dt)
            push!(pred_ts, t)
            upred .= p
            return nothing
        end
        prob = ODEProblem(recording_f!(evals), [1.0], (0.0, 1.0), [SEED])
        for make_alg in (
                nl -> ImplicitEuler(nlsolve = nl),
                nl -> TRBDF2(nlsolve = nl),
                nl -> FBDF(nlsolve = nl),
            )
            empty!(evals)
            empty!(pred_ts)
            sol = solve(
                prob, make_alg(NLNewton(predictor = pred!));
                dt = 0.1, adaptive = false
            )
            @test sol.retcode == SciMLBase.ReturnCode.Success
            nsteps = length(sol.t) - 1
            @test length(pred_ts) >= nsteps
            @test count(==(SEED), evals) >= nsteps
        end
    end

    @testset "all nonlinear solver backends" begin
        nlsolve_algs = [
            NLNewton,
            NLFunctional,
            NLAnderson,
            (; kwargs...) -> NonlinearSolveAlg(NewtonRaphson(); kwargs...),
        ]
        for make_nlsolve in nlsolve_algs
            evals = Float64[]
            pred = (uprev, p, t, dt) -> SEED
            prob = ODEProblem(recording_f(evals), 1.0, (0.0, 1.0))
            sol = solve(
                prob, ImplicitEuler(nlsolve = make_nlsolve(predictor = pred));
                dt = 0.1, adaptive = false
            )
            @test sol.retcode == SciMLBase.ReturnCode.Success
            @test count(==(SEED), evals) >= 10
        end

        for make_nlsolve in nlsolve_algs
            evals = Float64[]
            pred! = (upred, uprev, p, t, dt) -> (upred .= SEED; nothing)
            prob = ODEProblem(recording_f!(evals), [1.0], (0.0, 1.0))
            sol = solve(
                prob, ImplicitEuler(nlsolve = make_nlsolve(predictor = pred!));
                dt = 0.1, adaptive = false
            )
            @test sol.retcode == SciMLBase.ReturnCode.Success
            @test count(==(SEED), evals) >= 10
        end
    end

    @testset "multi-stage methods seed every stage" begin
        # KenCarp4 is an ESDIRK (first stage explicit): the callable must fire
        # once per *implicit* stage, at distinct stage times.
        evals = Float64[]
        pred_ts = Float64[]
        pred = (uprev, p, t, dt) -> (push!(pred_ts, t); SEED)
        prob = ODEProblem(recording_f(evals), 1.0, (0.0, 1.0))
        sol = solve(
            prob, KenCarp4(nlsolve = NLNewton(predictor = pred));
            dt = 0.1, adaptive = false
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test length(pred_ts) >= 2 * 10
        @test allunique(pred_ts)
        @test count(==(SEED), evals) >= 2 * 10
    end

    @testset "a good predictor reduces nonlinear iterations" begin
        # u' = u with ImplicitEuler: the converged stage value is uprev/(1 - dt),
        # so this predictor lands exactly on the solution of each stage solve.
        pred = (uprev, p, t, dt) -> uprev / (1 - dt)
        prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0))
        sol_pred = solve(
            prob, ImplicitEuler(nlsolve = NLNewton(predictor = pred));
            dt = 0.1, adaptive = false
        )
        sol_def = solve(prob, ImplicitEuler(); dt = 0.1, adaptive = false)
        @test sol_pred.retcode == SciMLBase.ReturnCode.Success
        @test sol_pred.u[end] ≈ sol_def.u[end] atol = 1.0e-8
        @test sol_pred.stats.nnonliniter < sol_def.stats.nnonliniter
        @test sol_pred.stats.nnonliniter <= 2 * 10
    end

    @testset "callable overrides the algorithm's Predictor enum" begin
        # With a callable set, `nlsolve!` overwrites whatever the kernel seeded —
        # and the kernel skips the enum machinery entirely (no `addsteps!` for
        # `MaxOrder`). The seed must still land for each of them.
        for enum_predictor in (Predictor.Linear, Predictor.MaxOrder)
            evals = Float64[]
            pred = (uprev, p, t, dt) -> SEED
            prob = ODEProblem(recording_f(evals), 1.0, (0.0, 1.0))
            sol = solve(
                prob,
                ImplicitEuler(
                    predictor = enum_predictor,
                    nlsolve = NLNewton(predictor = pred)
                );
                dt = 0.1, adaptive = false
            )
            @test sol.retcode == SciMLBase.ReturnCode.Success
            @test count(==(SEED), evals) >= 10
        end
    end

    @testset "rejects misuse" begin
        # The `Predictor` enum belongs to the ODE algorithm, not the nlsolve alg.
        @test_throws ArgumentError NLNewton(predictor = Predictor.Linear)
        @test_throws ArgumentError NLFunctional(predictor = Predictor.Trivial)
        @test_throws ArgumentError NLAnderson(predictor = Predictor.MaxOrder)
        @test_throws ArgumentError NonlinearSolveAlg(
            NewtonRaphson(); predictor = Predictor.Trivial
        )
        # Homotopy solvers start from the λ = 0 anchor, not an initial guess.
        @test_throws ArgumentError HomotopyNonlinearSolveAlg(
            predictor = (uprev, p, t, dt) -> uprev
        )
        @test HomotopyNonlinearSolveAlg() isa Any
    end

    @testset "DAE stage solves" begin
        # du[1] = u[2], u[1] + u[2] = exp(-t): index-1 semi-explicit DAE.
        evals = Vector{Float64}[]
        function dae!(out, du, u, p, t)
            eltype(u) === Float64 && push!(evals, copy(u))
            out[1] = du[1] - u[2]
            out[2] = u[1] + u[2] - exp(-t)
            return nothing
        end
        u0 = [0.0, 1.0]
        du0 = [1.0, 0.0]
        prob = DAEProblem(dae!, du0, u0, (0.0, 0.5); differential_vars = [true, false])
        pred! = (upred, uprev, p, t, dt) -> (upred .= SEED; nothing)
        sol = solve(
            prob, DImplicitEuler(nlsolve = NLNewton(predictor = pred!));
            dt = 0.1, adaptive = false
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test any(==([SEED, SEED]), evals)
    end
end

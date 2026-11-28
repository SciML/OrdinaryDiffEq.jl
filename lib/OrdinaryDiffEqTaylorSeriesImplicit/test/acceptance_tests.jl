module TaylorAcceptanceTests

using Test, LinearAlgebra, LinearSolve, OrdinaryDiffEqTaylorSeriesImplicit
using SciMLBase: SciMLBase, ReturnCode, successful_retcode

const ODTI = OrdinaryDiffEqTaylorSeriesImplicit

# The [2/3] Padé equation for this one step is (8/5)u - 2/15 = 0.
# At the initial guess the correction to u cancels while the stages are still wrong.
function cancelling_stages!(du, u, p, t)
    @. du = 4 * (t - 1) * u + 3t^2 - 2t
    return nothing
end

decay!(du, u, p, t) = (du .= -u; nothing)

struct FaultSolver <: LinearSolve.AbstractSolveFunction
    target::Base.RefValue{Any}
    nonfinite::Bool
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::FaultSolver; kwargs...)
    if cache === alg.target[]
        u = fill(alg.nonfinite ? NaN : 0.0, size(cache.u))
        code = alg.nonfinite ? ReturnCode.Success : ReturnCode.Failure
    else
        u = cache.A \ cache.b
        code = ReturnCode.Success
    end
    return SciMLBase.build_linear_solution(alg, u, nothing, cache; retcode = code)
end

@testset "Implicit Taylor independent acceptance" begin
    prob = ODEProblem(cancelling_stages!, [0.0], (0.0, 1.0))
    alg = TaylorRadauLA5(stages = :derivative, maxiters = 100)
    for adaptive in (false, true)
        it = init(prob, alg; dt = 1.0, adaptive, abstol = 1.0e-8, reltol = 1.0e-8)
        cache = ODTI.newton_state(it.cache)
        converged = ODTI.implicit_taylor_step!(it, cache, it.t, it.dt, cache.uprev, cache.u)
        if adaptive && !converged
            # An adaptive attempt may reject this long initial step.
            @test it.force_stepfail
            @test !cache.stages_valid
        else
            @test converged
            @test cache.u[1] ≈ 1 / 12 atol = 1.0e-8
        end
    end
    sol = solve(prob, alg; dt = 1.0, abstol = 1.0e-8, reltol = 1.0e-8)
    @test successful_retcode(sol)
    # ODE solution at t=1, from a high-precision power series of the affine RHS.
    @test sol.u[end][1] ≈ 0.08058583246643535 atol = 1.0e-6

    # The fixed-step Padé denominator is invertible at hλ = 5 even though
    # the adaptive error filter has a pole there.
    growth!(du, u, p, t) = (du .= 5 .* u; nothing)
    sol = solve(ODEProblem(growth!, [1.0], (0.0, 1.0)), alg; dt = 1.0, adaptive = false)
    @test successful_retcode(sol)
    @test sol.u[end][1] ≈ -51 / 4 atol = 1.0e-9
end

@testset "Failed linear solves reject the candidate" begin
    prob = ODEProblem(decay!, [1.0, 2.0], (0.0, 1.0))
    for (stages, location) in ((:reduced, :newton), (:reduced, :filter), (:derivative, :real), (:derivative, :complex)),
            nonfinite in (false, true)
        fault = Ref{Any}(nothing)
        solver = FaultSolver(fault, nonfinite)
        it = init(
            prob, TaylorRadauLA5(; stages, linsolve = solver);
            dt = 0.1, abstol = 1.0e-10, reltol = 1.0e-8
        )
        c = ODTI.newton_state(it.cache)
        fault[] = if location === :newton
            c.linsolve
        elseif location === :filter
            c.linsolve_filter
        elseif location === :real
            first(c.stage_factors.lr)
        else
            first(c.stage_factors.lc)
        end
        @test !ODTI.implicit_taylor_step!(it, c, it.t, it.dt, c.uprev, c.u)
        @test it.force_stepfail
        @test !(c.lu_valid || c.lu_stages_valid || c.filter_valid || c.stages_valid)
        @test c.u == c.uprev
        fault[] = nothing
        SciMLBase.reinit!(it)
        sol = solve!(it)
        @test successful_retcode(sol)
        @test sol.u[end] ≈ exp(-1) .* prob.u0 rtol = 1.0e-5
    end

    # A zero first pivot leaves the second pivot uninitialized.
    growth!(du, u, p, t) = (du .= u; nothing)
    sol = solve(ODEProblem(growth!, [1.0, 2.0], (0.0, 1.0)), ImplicitTaylor(); dt = 1.0, adaptive = false)
    @test sol.retcode == ReturnCode.ConvergenceFailure
end

@testset "Derivative factors preserve Krylov right-hand sides" begin
    it = init(
        ODEProblem(decay!, [1.0, 2.0], (0.0, 1.0)),
        TaylorRadauLA5(stages = :derivative, linsolve = KrylovJL_GMRES())
    )
    S = ODTI.newton_state(it.cache).stage_factors
    for (matrices, solvers, buffers) in ((S.Ar, S.lr, S.wr), (S.Ac, S.lc, S.wc))
        for (A, solver, w) in zip(matrices, solvers, buffers)
            A .= [2.0 0.25; -0.5 3.0]
            eltype(A) <: Complex && (A[1, 1] += 0.7im)
            original = copy(A)
            for k in 1:3
                b = eltype(w)[k, -0.3k]
                copyto!(w, b)
                ODTI.solve_stage_factor!(it, solver, A, w, k == 1)
                @test w ≈ original \ b rtol = 1.0e-8
            end
        end
    end
end

end

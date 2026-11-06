using OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg, tracks_J_age
using NonlinearSolve: NewtonRaphson
using ODEProblemLibrary
using ADTypes, SciMLBase
using Test

nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))

@testset "a step the error test rejected does not reuse its W" begin
    # The controller shrinks `dt` by ~15% after a rejected step, which leaves `γΔt` inside
    # `new_W_dt_cutoff` (1//5) — so the retry reused a `W` whose chord iteration contracts at
    # ~0.25 per iteration and cannot converge inside `maxiters`. Every such stage was charged
    # to the step controller as a nonlinear convergence failure. `do_newJW` has always
    # reassembled `W` there; `NonlinearSolveAlg` did not, and reported 46-135 failures where
    # `NLNewton` reports 0-1.
    for (pname, prob) in (
            ("brusselator1d", ODEProblemLibrary.prob_ode_brusselator_1d),
            ("hires", ODEProblemLibrary.prob_ode_hires),
            ("pollution", ODEProblemLibrary.prob_ode_pollution),
        )
        ref = solve(prob, TRBDF2(); reltol = 1.0e-9, abstol = 1.0e-12)
        sol = solve(prob, TRBDF2(nlsolve = nsa); reltol = 1.0e-9, abstol = 1.0e-12)
        @testset "$pname" begin
            @test SciMLBase.successful_retcode(sol)
            @test sol.stats.nnonlinconvfail <= ref.stats.nnonlinconvfail + 2
            @test sol.stats.njacs <= 4 * max(ref.stats.njacs, 1)
            @test sol.stats.nf < 1.1 * ref.stats.nf
        end
    end
end

@testset "only caches that record the Jacobian's age take part in the retry" begin
    # The retry re-enters `nlsolve!` until `isJcurrent`, so a cache with no reusable `W` to
    # refresh — and hence no `J_t` to advance — must stay out of it or the solve would hang.
    lin(u, p, t) = [-1.0e4 * u[1] + u[2], -u[2]]
    function lin!(du, u, p, t)
        du[1] = -1.0e4 * u[1] + u[2]
        du[2] = -u[2]
        return nothing
    end
    for (label, f, tracks) in (
            ("in-place, reusable W", lin!, true),
            ("out-of-place, no W", lin, false),
        )
        prob = ODEProblem(f, [1.0, 1.0], (0.0, 1.0))
        integ = init(prob, TRBDF2(nlsolve = nsa); reltol = 1.0e-8, abstol = 1.0e-10)
        @testset "$label" begin
            @test tracks_J_age(integ.cache.nlsolver) == tracks
            solve!(integ)
            @test SciMLBase.successful_retcode(integ.sol)
        end
    end
end

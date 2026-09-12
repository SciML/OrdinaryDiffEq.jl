using OrdinaryDiffEqNonlinearSolve: NLNewton, NonlinearSolveAlg
using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, SciMLBase
using Test

@testset "stage predictor on the nonlinear solver" begin
    calls = Ref(0)
    seed_oop = (uprev, p, t, dt) -> (calls[] += 1; uprev)
    seed_iip = (upred, uprev, p, t, dt) -> (calls[] += 1; upred .= uprev)
    prob_oop = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    prob_iip = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 2.0], (0.0, 1.0))

    @testset "$(nameof(M)) $(nameof(N)) $(iip ? "iip" : "oop")" for M in
            (ImplicitEuler, TRBDF2, FBDF), N in (NLNewton, NonlinearSolveAlg), iip in (false, true)

        prob = iip ? prob_iip : prob_oop
        seed = iip ? seed_iip : seed_oop
        calls[] = 0
        guessed = solve(prob, M(nlsolve = N(predictor = seed)); abstol = 1.0e-10, reltol = 1.0e-10)
        @test calls[] > 0
        @test SciMLBase.successful_retcode(guessed)
        plain = solve(prob, M(nlsolve = N()); abstol = 1.0e-10, reltol = 1.0e-10)
        @test guessed.u[end] ≈ plain.u[end] rtol = 1.0e-6
    end
end

using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using ADTypes, LinearAlgebra, SciMLBase
using Test

# Van der Pol μ=1e5: the smoothed (W⁻¹-filtered) SDIRK error estimate and the raw
# embedded estimate differ by orders of magnitude in accepted step counts here, which
# makes it a sharp detector for whether the smoothing is actually active. Before the
# fix the smoothing was `isnewton`-gated, so NonlinearSolveAlg silently error-controlled
# on the raw estimate: `smooth_est` had no effect and step counts diverged ~40x from
# NLNewton on the same `TRBDF2()` call.
function vdp!(du, u, p, t)
    du[1] = u[2]
    du[2] = p[1] * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end
prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [1.0e5])
nsa() = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))

@testset "smooth_est is active under NonlinearSolveAlg W-reuse" begin
    for ALG in (TRBDF2, KenCarp4)
        s_smooth = solve(
            prob, ALG(nlsolve = nsa()); reltol = 1.0e-8, abstol = 1.0e-11
        )
        s_raw = solve(
            prob, ALG(smooth_est = false, nlsolve = nsa());
            reltol = 1.0e-8, abstol = 1.0e-11
        )
        @test SciMLBase.successful_retcode(s_smooth)
        @test SciMLBase.successful_retcode(s_raw)
        # smoothing must change the error control (pre-fix these were identical)
        @test s_smooth.stats.naccept < s_raw.stats.naccept / 2

        # and NSA must now track NLNewton under the same estimator settings
        s_nln = solve(prob, ALG(); reltol = 1.0e-8, abstol = 1.0e-11)
        @test s_smooth.stats.naccept < 3 * s_nln.stats.naccept
        s_nln_raw = solve(
            prob, ALG(smooth_est = false); reltol = 1.0e-8, abstol = 1.0e-11
        )
        @test 0.5 < s_raw.stats.naccept / s_nln_raw.stats.naccept < 2
    end
end

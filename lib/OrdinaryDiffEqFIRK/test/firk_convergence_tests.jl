using OrdinaryDiffEqFIRK, OrdinaryDiffEqCore, Test

# Stale Newton contraction estimate `η` made the stopping test
# `η * ndw < κ` pass on the first iteration regardless of `reltol`, so the solution
# error stopped responding to tolerance, as `RadauIIA9` returned 3.871e-7 at 1e-6,
# 1e-8 and 1e-10 alike.

# u_i(t) = exp(-t) solves this exactly. The coupling switching on at t = 0.2 shifts
# the Newton contraction mid-integration, which is what a stale `η` fails to notice.
function _tol_tracking_f!(du, u, p, t)
    e = exp(-t)
    for i in eachindex(u)
        du[i] = -10i * (u[i] - e) - e
        if t > 0.2 && i < length(u)
            du[i] += 3 * (u[i + 1] - e)
        end
    end
    return nothing
end

_tol_tracking_f(u, p, t) = (du = similar(u); _tol_tracking_f!(du, u, p, t); du)

const TOL_TRACKING_ALGS = (RadauIIA3(), RadauIIA5(), RadauIIA9(), AdaptiveRadau())
const TOL_TRACKING_TOLS = (1.0e-6, 1.0e-8, 1.0e-10)

"""Endpoint error of `prob` under `alg` at relative tolerance `tol`."""
function tol_tracking_error(prob, alg, tol)
    sol = solve(
        prob, alg; reltol = tol, abstol = tol / 100, dt = 0.01, tstops = [0.2]
    )
    @test sol.retcode == ReturnCode.Success
    return maximum(abs.(sol.u[end] .- exp(-0.5)))
end

@testset "FIRK tolerance tracking (#4783)" begin
    for (tag, f) in (("iip", _tol_tracking_f!), ("oop", _tol_tracking_f))
        prob = ODEProblem(f, ones(8), (0.0, 0.5))
        @testset "$tag" begin
            for alg in TOL_TRACKING_ALGS
                errs = [tol_tracking_error(prob, alg, tol) for tol in TOL_TRACKING_TOLS]
                @test errs[end] < errs[1] / 10
                @test errs[end] < 1.0e-8
            end
        end
    end
end

using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF: NordsieckBDF
using ADTypes: AutoForwardDiff
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

# The first Newton iteration accepts on the rate `prev_θ` of an earlier solve. A
# rate measured with an old W must not carry over to a new one: here a stale small
# rate let every step stop after one iteration, leaving an error that flips sign
# each step and stalls the step size (SciML/OrdinaryDiffEq.jl#4983).
relaxed_sqrt(x, t) = abs(x) < t ? 1 / 4 * (x / sqrt(t)) * (5 - (x / t)^2) :
    sign(x) * sqrt(abs(x))
pit_level(u, t) = -0.0164 + (u[1] + u[2] - 1.0 * t) / 119.0
function pit!(du, u, p, t)
    h = pit_level(u, t)
    du[1] = 36.5 * relaxed_sqrt(0.0004 - h, 1.0e-4)
    du[2] = 36.5 * relaxed_sqrt(0.0 - h, 1.0e-4)
    for k in 3:12
        du[k] = (5.0 * t - u[k]) / 0.05
    end
    return nothing
end

@testset "prev_θ restarts on a new W" begin
    prob = ODEProblem(pit!, zeros(12), (0.0, 3600.0))
    alg = NordsieckBDF(; nlsolve = NLNewton(), autodiff = AutoForwardDiff())
    # RMS norm of a system with 2000 states, the others carrying no error
    rms = (u, t) -> u isa Number ? abs(u) : sqrt(sum(abs2, u) / 2000)
    sol = solve(
        prob, alg; abstol = 0.01, reltol = 0.001, internalnorm = rms,
        save_everystep = false
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.naccept < 100  # 2793 with the stale rate
    @test abs(pit_level(sol.u[end], 3600.0) + 0.05e-3) < 1.0e-4  # converged level
end

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    # BackTracking should cost at most marginally more than plain Newton
    @test sol2.stats.nf <= 1.001 * sol1.stats.nf
end

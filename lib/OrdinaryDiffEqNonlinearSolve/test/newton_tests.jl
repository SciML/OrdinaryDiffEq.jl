using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF: NordsieckBDF
using ADTypes: AutoForwardDiff
using SciMLOperators: MatrixOperator
using LinearAlgebra: Tridiagonal
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

# The first Newton iteration accepts on the rate `prev_θ` of an earlier solve. Once
# the pit level reaches the steep part of the square root, the Jacobian of the first
# step is far off. With a small rate remembered, every step still stopped after one
# iteration, so nothing noticed: the level overshot each step and the step size
# stalled (SciML/OrdinaryDiffEq.jl#4983). Restarting the rate on a new W forces a
# second iteration, which sees the divergence and refreshes the Jacobian.
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
    # The level settles where the two outflows carry off the inflow of 1. The stale
    # rate left it centimeters away, with the flows reversed.
    du = similar(sol.u[end])
    pit!(du, sol.u[end], nothing, 3600.0)
    @test du[1] + du[2] ≈ 1 rtol = 1.0e-2
end

@testset "a linear function keeps its rate" begin
    # Diffusion as a linear operator, treated implicitly, plus an explicit reaction.
    n = 40
    A = Matrix(Tridiagonal(fill(1.0, n - 1), fill(-2.0, n), fill(1.0, n - 1)))
    A .*= (n + 1)^2 / 10
    reaction!(du, u, p, t) = (@. du = 0.5 * cos(u); nothing)
    u0 = [sin(π * i / (n + 1)) for i in 1:n]
    prob = SplitODEProblem(MatrixOperator(A), reaction!, u0, (0.0, 2.0))
    sol = solve(prob, KenCarp4(); abstol = 1.0e-6, reltol = 1.0e-6)
    @test SciMLBase.successful_retcode(sol)
    # W follows the step size, but the Jacobian is exact, so each of the 5 implicit
    # stages of a step needs one iteration once the rate is known. Restarting the
    # rate on every new W takes two.
    nstages = 5 * (sol.stats.naccept + sol.stats.nreject)
    @test sol.stats.nnonliniter < 1.5 * nstages
end

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    # `stats.nf` counts the residual evaluations of the Newton iterations, not the ones
    # the line search adds, so this compares iteration counts and not work. BackTracking
    # takes the full Newton step in over 99.9% of the iterations here. The few it
    # shortens change the step sequence, which moves `nf` by up to 1.5% either way.
    # Shortening steps as a rule costs far more: `relax = 1 // 10` adds 70% to 100%.
    @test sol2.stats.nf <= 1.1 * sol1.stats.nf
end

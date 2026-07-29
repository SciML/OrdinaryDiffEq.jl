using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using ADTypes, LinearAlgebra, SciMLBase
using Test

function rober!(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    return nothing
end
const JAC_CALLS = Ref(0)
function rober_jac!(J, u, p, t)
    JAC_CALLS[] += 1
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    J[1, 1] = -k₁
    J[1, 2] = k₃ * y₃
    J[1, 3] = k₃ * y₂
    J[2, 1] = k₁
    J[2, 2] = -2k₂ * y₂ - k₃ * y₃
    J[2, 3] = -k₃ * y₂
    J[3, 1] = 0
    J[3, 2] = 2k₂ * y₂
    J[3, 3] = 0
    return nothing
end
f = ODEFunction(rober!; jac = rober_jac!)
prob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e5), [0.04, 3.0e7, 1.0e4])
nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))
refsol = solve(prob, FBDF(); reltol = 1.0e-12, abstol = 1.0e-14)

@testset "dt-only W updates reuse the stored Jacobian" begin
    for alg in (FBDF(nlsolve = nsa), TRBDF2(nlsolve = nsa))
        JAC_CALLS[] = 0
        sol = solve(prob, alg; reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        @test norm(sol.u[end] .- refsol.u[end]) / norm(refsol.u[end]) < 1.0e-4
        # W is reassembled on dt changes throughout (stats.nw ~ O(100) here), but a
        # fresh Jacobian evaluation is only needed on first use and after divergences.
        # Before the J/W split every W update called the user jac, so this tracked
        # stats.nw one-to-one.
        @test JAC_CALLS[] < sol.stats.nw / 3
        @test JAC_CALLS[] >= 1
    end
end

@testset "stale-W Newton directions are rescaled" begin
    # `W` is reused while `γΔt` drifts within `new_W_dt_cutoff`, which leaves the Newton
    # direction it produces scaled wrong; `NLNewton` corrects for that. Reproducing it
    # needs `dt` to range over many orders of magnitude, hence the long tspan, and an
    # AD Jacobian rather than the analytic one used above.
    staleprob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e11), [0.04, 3.0e7, 1.0e4])
    tight = solve(staleprob, Rodas5P(); reltol = 1.0e-13, abstol = 1.0e-15)
    err(sol) = maximum(abs.(sol.u[end] .- tight.u[end]))

    ref = solve(staleprob, KenCarp4(); reltol = 1.0e-9, abstol = 1.0e-12)
    sol = solve(staleprob, KenCarp4(nlsolve = nsa); reltol = 1.0e-9, abstol = 1.0e-12)
    @test SciMLBase.successful_retcode(sol)
    # Uncorrected directions cost ~1.5x NLNewton's steps for ~9x its error here.
    @test err(sol) < 3 * err(ref)
    @test sol.stats.naccept < 1.25 * ref.stats.naccept
end

using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
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

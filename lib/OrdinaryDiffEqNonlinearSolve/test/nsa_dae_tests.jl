using OrdinaryDiffEqBDF
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using ADTypes, LinearAlgebra, SciMLBase
using Test

# Robertson in fully implicit residual form, third equation algebraic. The stage residual of
# a `DAEProblem` goes through `daenlf`/`oopdaenlf`, which are the only callers of the
# `DAEFunction` methods of `_compute_rhs!`/`_compute_rhs` — a signature the adapters used to
# disagree with, so every `NonlinearSolveAlg` call on a true DAE raised a `MethodError`.
function rober_dae!(out, du, u, p, t)
    k₁, k₂, k₃ = p
    out[1] = -k₁ * u[1] + k₃ * u[2] * u[3] - du[1]
    out[2] = k₁ * u[1] - k₂ * u[2]^2 - k₃ * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
    return nothing
end

function rober_dae(du, u, p, t)
    k₁, k₂, k₃ = p
    return [
        -k₁ * u[1] + k₃ * u[2] * u[3] - du[1],
        k₁ * u[1] - k₂ * u[2]^2 - k₃ * u[2] * u[3] - du[2],
        u[1] + u[2] + u[3] - 1.0,
    ]
end

# Equivalent singular-mass-matrix ODE, for a reference independent of the DAE solvers.
function rober_mm!(du, u, p, t)
    k₁, k₂, k₃ = p
    du[1] = -k₁ * u[1] + k₃ * u[2] * u[3]
    du[2] = k₁ * u[1] - k₂ * u[2]^2 - k₃ * u[2] * u[3]
    du[3] = u[1] + u[2] + u[3] - 1.0
    return nothing
end

const P = [0.04, 3.0e7, 1.0e4]
const U0 = [1.0, 0.0, 0.0]
const DU0 = [-0.04, 0.04, 0.0]
const TSPAN = (0.0, 1.0e4)
const DVARS = [true, true, false]

prob_iip = DAEProblem(rober_dae!, DU0, U0, TSPAN, P; differential_vars = DVARS)
prob_oop = DAEProblem(
    DAEFunction{false}(rober_dae), DU0, U0, TSPAN, P; differential_vars = DVARS
)
prob_mm = ODEProblem(
    ODEFunction(rober_mm!; mass_matrix = [1.0 0 0; 0 1.0 0; 0 0 0]), U0, TSPAN, P
)
uref = solve(prob_mm, Rodas5P(); reltol = 1.0e-12, abstol = 1.0e-12).u[end]

nsa(ad) = NonlinearSolveAlg(NewtonRaphson(; autodiff = ad))

@testset "NonlinearSolveAlg solves a true DAEProblem" begin
    @testset "$pname / $aname / $adname" for (pname, prob) in
            (("in-place", prob_iip), ("out-of-place", prob_oop)),
            (aname, ALG) in (("DFBDF", DFBDF), ("DImplicitEuler", DImplicitEuler)),
            (adname, ad) in
            (("AutoFiniteDiff", AutoFiniteDiff()), ("AutoForwardDiff", AutoForwardDiff()))

        sol = solve(prob, ALG(nlsolve = nsa(ad)); reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(sol)

        ref = solve(prob, ALG(); reltol = 1.0e-8, abstol = 1.0e-8)
        @test SciMLBase.successful_retcode(ref)

        # Right answer, not merely a successful retcode: within the accuracy `NLNewton`
        # itself reaches on this problem at this tolerance.
        tol = 4 * maximum(abs, (ref.u[end] .- uref) ./ max.(abs.(uref), 1.0e-8)) + 1.0e-6
        @test maximum(abs, (sol.u[end] .- uref) ./ max.(abs.(uref), 1.0e-8)) < tol

        # And it must get there by taking a comparable number of steps: a stage that never
        # converges still integrates the problem, just with `dt` driven to nothing.
        @test 0.25 < sol.stats.naccept / ref.stats.naccept < 4
    end
end

using OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg, NLNewton
using NonlinearSolve: NewtonRaphson
using ADTypes, SciMLBase
using Test

const FCALLS = Ref(0)
const JCALLS = Ref(0)

function rober!(du, u, p, t)
    FCALLS[] += 1
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    du[3] = k₂ * y₂^2
    return nothing
end

function rober_jac!(J, u, p, t)
    JCALLS[] += 1
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

prob = ODEProblem(
    ODEFunction(rober!; jac = rober_jac!), [1.0, 0.0, 0.0],
    (0.0, 1.0e5), [0.04, 3.0e7, 1.0e4]
)
nsa = NonlinearSolveAlg(NewtonRaphson(; autodiff = AutoForwardDiff()))

@testset "reported work matches the work actually performed" begin
    # `stats` is what every work-precision comparison reads, so it has to count the same
    # quantity for both nonlinear solvers. `NonlinearSolveAlg` used to forward the inner
    # cache's counters verbatim: that missed the residual evaluation `reinit!` performs
    # before zeroing them, and counted `WReuseJac` copies as Jacobian evaluations.
    for A in (KenCarp4, TRBDF2, FBDF), nl in (NLNewton(), nsa)
        FCALLS[] = 0
        JCALLS[] = 0
        sol = solve(prob, A(nlsolve = nl); reltol = 1.0e-8, abstol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        # Never claim more work than was done, and account for all but a handful of
        # calls (cache construction is outside the per-stage accounting).
        @test sol.stats.nf <= FCALLS[]
        @test sol.stats.nf >= 0.99 * FCALLS[]
        @test sol.stats.njacs == JCALLS[]
    end
end

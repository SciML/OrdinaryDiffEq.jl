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
const LINJAC = [-1.5 0.25; 0.0 -0.5]

function lin!(du, u, p, t)
    FCALLS[] += 1
    du[1] = LINJAC[1, 1] * u[1] + LINJAC[1, 2] * u[2]
    du[2] = LINJAC[2, 2] * u[2]
    return nothing
end

linprob = ODEProblem(
    ODEFunction(lin!; jac = (J, u, p, t) -> (J .= LINJAC; nothing)),
    [1.0, 1.0], (0.0, 20.0)
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

@testset "no residual evaluation per stage beyond the Newton iterations" begin
    # A `k`-iteration stage needs the residual at the `k` points its Newton steps are taken
    # from, but a `NonlinearSolve` cache evaluates it where each step *lands*, so driving it
    # naively costs `k + 1`. The Jacobians here are analytic, so `NLNewton` measures how many
    # `f` calls outside the Newton iterations a step legitimately needs.
    #
    # `linprob` is the case where deferring buys nothing and so has to cost nothing: its
    # stage equation is solved exactly by one Newton step. Only one-implicit-stage methods
    # are compared on it, because the excess below is per step while `NonlinearSolveAlg`
    # spends one residual per *stage* re-initializing the inner cache.
    for (p, algs) in ((prob, (KenCarp4, TRBDF2, FBDF)), (linprob, (FBDF, QNDF))), A in algs
        excess = map((NLNewton(), nsa)) do nl
            FCALLS[] = 0
            sol = solve(p, A(nlsolve = nl); reltol = 1.0e-8, abstol = 1.0e-10)
            @test SciMLBase.successful_retcode(sol)
            (FCALLS[] - sol.stats.nnonliniter) / sol.stats.naccept
        end
        @test excess[2] <= excess[1] + 0.25
    end
end

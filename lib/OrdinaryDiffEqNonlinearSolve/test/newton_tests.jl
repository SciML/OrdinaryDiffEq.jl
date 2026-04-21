using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using LinearSolve
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end

# Canonical linsolve location is on the nlsolve object. Setting it on the nlsolve
# should be equivalent to setting it via the legacy alg-level kwarg, and when
# both are set the nlsolve-level linsolve wins.
@testset "NLNewton linsolve field" begin
    prob = prob_ode_orego
    ref = solve(prob, Trapezoid(), reltol = 1.0e-10, abstol = 1.0e-10)

    via_alg = solve(
        prob, Trapezoid(linsolve = LUFactorization()),
        reltol = 1.0e-10, abstol = 1.0e-10
    )
    via_nlsolve = solve(
        prob, Trapezoid(nlsolve = NLNewton(linsolve = LUFactorization())),
        reltol = 1.0e-10, abstol = 1.0e-10
    )
    @test via_alg.retcode == SciMLBase.ReturnCode.Success
    @test via_nlsolve.retcode == SciMLBase.ReturnCode.Success
    @test via_alg.u ≈ ref.u
    @test via_nlsolve.u ≈ ref.u

    # nlsolve-level linsolve wins over the alg-level fallback
    both = solve(
        prob,
        Trapezoid(
            linsolve = KLUFactorization(),
            nlsolve = NLNewton(linsolve = LUFactorization())
        ),
        reltol = 1.0e-10, abstol = 1.0e-10
    )
    @test both.retcode == SciMLBase.ReturnCode.Success
    @test both.u ≈ ref.u
end

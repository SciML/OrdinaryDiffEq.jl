using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
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

# Forged non-negative FD slope at α=0: BackTracking throws, linesearch_step salvages.
@testset "non-negative FD slope is salvaged for BackTracking" begin
    ϕ(α) = abs2(one(α) - α)
    dϕ(α) = oftype(α, 2) * (α - one(α))
    ϕdϕ(α) = iszero(α) ? (ϕ(α), one(α)) : (ϕ(α), dϕ(α))

    ϕ0, dϕ0_fd = ϕdϕ(zero(1.0))
    @test dϕ0_fd > 0
    @test_throws LineSearches.LineSearchException begin
        BackTracking()(ϕ, dϕ, ϕdϕ, one(ϕ0), ϕ0, dϕ0_fd)
    end

    α = OrdinaryDiffEqNonlinearSolve.linesearch_step(
        BackTracking(), ϕ, dϕ, ϕdϕ, one(ϕ0)
    )
    @test α > 0
    @test ϕ(α) < ϕ0
end

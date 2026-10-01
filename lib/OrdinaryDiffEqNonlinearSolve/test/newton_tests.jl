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

# Forged non-negative FD slope at α=0 with ϕ(1)>ϕ0: Kelley salvage must backtrack.
@testset "non-negative FD slope is salvaged for BackTracking" begin
    ϕ(α) = abs2(one(α) - oftype(α, 4) * α)
    dϕ(α) = oftype(α, -8) * (one(α) - oftype(α, 4) * α)
    ϕdϕ(α) = iszero(α) ? (ϕ(α), one(α)) : (ϕ(α), dϕ(α))

    ϕ0 = ϕ(zero(1.0))
    @test ϕ(one(ϕ0)) > ϕ0

    α = OrdinaryDiffEqNonlinearSolve.linesearch_step(
        BackTracking(), ϕ, dϕ, ϕdϕ, one(ϕ0)
    )
    @test α < one(α)
    @test ϕ(α) < ϕ0
end

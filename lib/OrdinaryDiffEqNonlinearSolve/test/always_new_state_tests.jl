using OrdinaryDiffEqSDIRK: TRBDF2
using OrdinaryDiffEqNonlinearSolve: NLNewton
using SciMLBase: ODEProblem, solve, successful_retcode
using Test

@testset "Always-new Newton preserves the previous state" begin
    rhs(u, p, t) = -u .* u
    function rhs!(du, u, p, t)
        @. du = -u * u
        return nothing
    end
    alg = TRBDF2(nlsolve = NLNewton(always_new = true))
    inplace = solve(
        ODEProblem(rhs!, ones(2), (0.0, 1.0)), alg; dt = 0.01, adaptive = false
    )
    outofplace = solve(
        ODEProblem(rhs, ones(2), (0.0, 1.0)), alg; dt = 0.01, adaptive = false
    )
    @test successful_retcode(inplace)
    @test successful_retcode(outofplace)
    @test inplace.u[end] ≈ fill(0.5, 2) atol = 1.0e-5
    @test outofplace.u[end] ≈ fill(0.5, 2) atol = 1.0e-5
    @test outofplace.u[end] ≈ inplace.u[end] rtol = 1.0e-10
end

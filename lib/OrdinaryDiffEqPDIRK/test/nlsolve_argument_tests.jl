using OrdinaryDiffEqPDIRK
using OrdinaryDiffEqCore
using OrdinaryDiffEqNonlinearSolve: nlsolve!
using SciMLBase
using Test

# The out-of-place branch used to pass `γ * dt` where `nlsolve!` takes the
# integrator cache. Nothing complained: the `NLNewton` guard only rejected
# `nothing`, and `update_W!` happens not to read that argument on the
# out-of-place path. Both halves are covered here.
@testset "PDIRK44 out-of-place agrees with in-place" begin
    oop = ODEProblem((u, p, t) -> -u * (1 + 0.1u), 1.0, (0.0, 1.0))
    iip = ODEProblem((du, u, p, t) -> (du[1] = -u[1] * (1 + 0.1u[1]); nothing), [1.0], (0.0, 1.0))
    for threading in (false, true)
        a = solve(oop, PDIRK44(; threading); dt = 0.05, adaptive = false)
        b = solve(iip, PDIRK44(; threading); dt = 0.05, adaptive = false)
        @test SciMLBase.successful_retcode(a)
        @test a.u[end] ≈ b.u[end][1] rtol = 1.0e-10
    end
end

@testset "nlsolve! rejects a third argument that is not a cache" begin
    prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))
    integrator = init(prob, PDIRK44(threading = false); dt = 0.1, adaptive = false)
    nlsolver = first(integrator.cache.nlsolver)
    @test_throws ArgumentError nlsolve!(nlsolver, integrator, 0.5, false)
    @test_throws ArgumentError nlsolve!(nlsolver, integrator, nothing, false)
end

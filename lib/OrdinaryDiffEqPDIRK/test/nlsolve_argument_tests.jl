using OrdinaryDiffEqPDIRK
using SciMLBase
using Test

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

using OrdinaryDiffEqPDIRK
using SciMLBase
using Test

@testset "PDIRK44 exposes scratch buffers from its solver vector" begin
    prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 1.0], (0.0, 1.0))
    integrator = init(prob, PDIRK44(); dt = 0.1, adaptive = false)
    step!(integrator)

    @test integrator.cache.nlsolver isa AbstractVector

    tmp_cache = SciMLBase.get_tmp_cache(integrator)
    @test tmp_cache isa NTuple{2, <:AbstractVector}
    @test all(buf -> axes(buf) == axes(prob.u0), tmp_cache)
end

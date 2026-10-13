using StochasticDiffEqHighOrder
using SciMLBase
using Test

# SRA() builds an SRACache with u_cache/du_cache/user_cache methods. Those must
# extend SciMLBase's functions, not HighOrder-local ones (otherwise the umbrella
# reexport clashes / leaves undefined exports).
@testset "SRACache extends SciMLBase cache hooks" begin
    @test StochasticDiffEqHighOrder.du_cache === SciMLBase.du_cache
    @test StochasticDiffEqHighOrder.u_cache === SciMLBase.u_cache
    @test StochasticDiffEqHighOrder.user_cache === SciMLBase.user_cache

    f!(du, u, p, t) = (du .= u)
    g!(du, u, p, t) = (du .= u)
    prob = SDEProblem(f!, g!, ones(2), (0.0, 1.0))
    integ = init(prob, SRA(); dt = 0.1)
    cache = integ.cache
    @test cache isa StochasticDiffEqHighOrder.SRACache

    @test !isempty(methods(SciMLBase.du_cache, (typeof(cache),)))
    @test !isempty(methods(SciMLBase.u_cache, (typeof(cache),)))
    @test !isempty(methods(SciMLBase.user_cache, (typeof(cache),)))

    du_bufs = SciMLBase.du_cache(cache)
    @test !isempty(du_bufs)
    user_bufs = SciMLBase.user_cache(cache)
    @test integ.u in user_bufs
end

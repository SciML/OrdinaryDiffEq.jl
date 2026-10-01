using StochasticDiffEqHighOrder
using SciMLBase
using Test

# SRA() builds an SRACache, which defines u_cache/du_cache/user_cache. Those
# methods must extend SciMLBase's functions (reached via DiffEqBase), not a
# HighOrder-local function, or resize! never sees them.
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

    du_bufs = SciMLBase.du_cache(cache)
    @test !isempty(du_bufs)
    user_bufs = SciMLBase.user_cache(cache)
    @test integ.u in user_bufs

    old_len = length(integ.u)
    resize!(integ, old_len + 1)
    @test length(integ.u) == old_len + 1
    for buf in du_bufs
        @test length(buf) == old_len + 1
    end
end

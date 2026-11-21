using OrdinaryDiffEqCore
using OrdinaryDiffEqCore: StatsDelta, stats_sink, charge_nf!, charge_nw!,
    charge_nsolve!, charge_njacs!, merge_stats_deltas!
using SciMLBase
using SciMLBase: DEStats
using Test

struct BufferedCache
    stats_delta::StatsDelta
end
struct PlainCache end
mutable struct FakeIntegrator
    stats::DEStats
end
SciMLBase.has_stats(::FakeIntegrator) = true

zeroed() = DEStats(0)

@testset "a cache without a buffer counts straight through" begin
    integrator = FakeIntegrator(zeroed())
    @test stats_sink(PlainCache()) === nothing
    @test stats_sink(nothing) === nothing
    charge_nf!(integrator, PlainCache(), 3)
    charge_nw!(integrator, PlainCache(), 1)
    @test integrator.stats.nf == 3
    @test integrator.stats.nw == 1
end

@testset "a cache with a buffer defers until drained" begin
    integrator = FakeIntegrator(zeroed())
    cache = BufferedCache(StatsDelta())
    @test stats_sink(cache) === cache.stats_delta

    charge_nf!(integrator, cache, 5)
    charge_nw!(integrator, cache, 2)
    charge_nsolve!(integrator, cache, 7)
    charge_njacs!(integrator, cache, 1)
    @test integrator.stats.nf == 0
    @test integrator.stats.nsolve == 0

    merge_stats_deltas!(integrator.stats, (cache,))
    @test integrator.stats.nf == 5
    @test integrator.stats.nw == 2
    @test integrator.stats.nsolve == 7
    @test integrator.stats.njacs == 1

    merge_stats_deltas!(integrator.stats, (cache,))
    @test integrator.stats.nf == 5
    @test integrator.stats.nsolve == 7
end

@testset "a StatsDelta is its own sink" begin
    integrator = FakeIntegrator(zeroed())
    delta = StatsDelta()
    @test stats_sink(delta) === delta
    charge_nf!(integrator, delta, 4)
    @test integrator.stats.nf == 0
    @test delta.nf == 4
end

@testset "workers accumulate privately and the parent merges after joining" begin
    integrator = FakeIntegrator(zeroed())
    caches = [BufferedCache(StatsDelta()) for _ in 1:64]
    Threads.@threads for cache in caches
        charge_nf!(integrator, cache, 10)
        charge_nsolve!(integrator, cache, 1)
    end
    @test integrator.stats.nf == 0
    merge_stats_deltas!(integrator.stats, caches)
    @test integrator.stats.nf == 640
    @test integrator.stats.nsolve == 64
end

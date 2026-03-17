using StochasticDiffEqLowOrder
using Test
using DiffEqBase, DiffEqNoiseProcess

@testset "StochasticDiffEqLowOrder" begin
    @testset "Module loads" begin
        @test isdefined(StochasticDiffEqLowOrder, :EM)
        @test isdefined(StochasticDiffEqLowOrder, :EulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEM)
    end
end

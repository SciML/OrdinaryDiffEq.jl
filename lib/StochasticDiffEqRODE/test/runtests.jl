using StochasticDiffEqRODE
using Test

@testset "StochasticDiffEqRODE" begin
    @testset "Module loads" begin
        @test isdefined(StochasticDiffEqRODE, :RandomEM)
        @test isdefined(StochasticDiffEqRODE, :RandomHeun)
        @test isdefined(StochasticDiffEqRODE, :RandomTamedEM)
        @test isdefined(StochasticDiffEqRODE, :BAOAB)
    end

    @testset "Algorithm construction" begin
        @test RandomEM() isa StochasticDiffEqRODEAlgorithm
        @test RandomHeun() isa StochasticDiffEqRODEAlgorithm
        @test RandomTamedEM() isa StochasticDiffEqRODEAlgorithm
        @test BAOAB() isa StochasticDiffEqAlgorithm
        @test BAOAB(gamma = 2.0, scale_noise = false).gamma == 2.0
        @test BAOAB(gamma = 2.0, scale_noise = false).scale_noise == false
    end
end

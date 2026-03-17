using StochasticDiffEqHighOrder
using Test
using DiffEqBase, DiffEqNoiseProcess

@testset "StochasticDiffEqHighOrder" begin
    @testset "Module loads" begin
        @test isdefined(StochasticDiffEqHighOrder, :SRI)
        @test isdefined(StochasticDiffEqHighOrder, :SRIW1)
        @test isdefined(StochasticDiffEqHighOrder, :SRIW2)
        @test isdefined(StochasticDiffEqHighOrder, :SOSRI)
        @test isdefined(StochasticDiffEqHighOrder, :SOSRI2)
        @test isdefined(StochasticDiffEqHighOrder, :SRA)
        @test isdefined(StochasticDiffEqHighOrder, :SRA1)
        @test isdefined(StochasticDiffEqHighOrder, :SRA2)
        @test isdefined(StochasticDiffEqHighOrder, :SRA3)
        @test isdefined(StochasticDiffEqHighOrder, :SOSRA)
        @test isdefined(StochasticDiffEqHighOrder, :SOSRA2)
    end
end

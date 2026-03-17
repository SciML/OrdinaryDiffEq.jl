using StochasticDiffEqMilstein
using Test
using DiffEqBase, DiffEqNoiseProcess

@testset "StochasticDiffEqMilstein" begin
    @testset "Module loads" begin
        @test isdefined(StochasticDiffEqMilstein, :RKMilGeneral)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_A)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_B)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_C)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_D)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_E)
        @test isdefined(StochasticDiffEqMilstein, :WangLi3SMil_F)
    end

    @testset "Algorithm constructors" begin
        @test RKMilGeneral() isa RKMilGeneral
        @test WangLi3SMil_A() isa WangLi3SMil_A
        @test WangLi3SMil_B() isa WangLi3SMil_B
        @test WangLi3SMil_C() isa WangLi3SMil_C
        @test WangLi3SMil_D() isa WangLi3SMil_D
        @test WangLi3SMil_E() isa WangLi3SMil_E
        @test WangLi3SMil_F() isa WangLi3SMil_F
    end
end

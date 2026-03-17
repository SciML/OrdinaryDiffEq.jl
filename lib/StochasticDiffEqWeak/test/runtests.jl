using StochasticDiffEqWeak
using Test

@testset "StochasticDiffEqWeak" begin
    @testset "Module loads" begin
        @test true
    end

    @testset "Algorithms are exported" begin
        @test DRI1 isa Type
        @test DRI1NM isa Type
        @test RI1 isa Type
        @test RI3 isa Type
        @test RI5 isa Type
        @test RI6 isa Type
        @test RDI1WM isa Type
        @test RDI2WM isa Type
        @test RDI3WM isa Type
        @test RDI4WM isa Type
        @test W2Ito1 isa Type
        @test RS1 isa Type
        @test RS2 isa Type
        @test PL1WM isa Type
        @test PL1WMA isa Type
        @test NON isa Type
        @test NON2 isa Type
        @test COM isa Type
        @test SIEA isa Type
        @test SIEB isa Type
        @test SMEA isa Type
        @test SMEB isa Type
        @test IRI1 isa Type
    end

    @testset "Algorithm constructors" begin
        @test DRI1() isa StochasticDiffEqWeak.StochasticDiffEqAdaptiveAlgorithm
        @test RDI1WM() isa StochasticDiffEqWeak.StochasticDiffEqAlgorithm
        @test RS1() isa StochasticDiffEqWeak.StochasticDiffEqAlgorithm
        @test NON() isa StochasticDiffEqWeak.StochasticDiffEqAlgorithm
        @test IRI1() isa StochasticDiffEqWeak.StochasticDiffEqNewtonAdaptiveAlgorithm
    end
end

using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqWeak
        using Test

        @test DRI1() isa StochasticDiffEqAdaptiveAlgorithm
        @test DRI1NM() isa StochasticDiffEqAdaptiveAlgorithm
        @test RI1() isa StochasticDiffEqAdaptiveAlgorithm
        @test RI3() isa StochasticDiffEqAdaptiveAlgorithm
        @test RI5() isa StochasticDiffEqAdaptiveAlgorithm
        @test RI6() isa StochasticDiffEqAdaptiveAlgorithm
        @test RDI1WM() isa StochasticDiffEqAlgorithm
        @test RDI2WM() isa StochasticDiffEqAdaptiveAlgorithm
        @test RDI3WM() isa StochasticDiffEqAdaptiveAlgorithm
        @test RDI4WM() isa StochasticDiffEqAdaptiveAlgorithm
        @test W2Ito1() isa StochasticDiffEqAdaptiveAlgorithm
        @test RS1() isa StochasticDiffEqAlgorithm
        @test RS2() isa StochasticDiffEqAlgorithm
        @test PL1WM() isa StochasticDiffEqAlgorithm
        @test PL1WMA() isa StochasticDiffEqAlgorithm
        @test NON() isa StochasticDiffEqAlgorithm
        @test NON2() isa StochasticDiffEqAlgorithm
        @test COM() isa StochasticDiffEqAlgorithm
        @test SIEA() isa StochasticDiffEqAlgorithm
        @test SIEB() isa StochasticDiffEqAlgorithm
        @test SMEA() isa StochasticDiffEqAlgorithm
        @test SMEB() isa StochasticDiffEqAlgorithm
        @test IRI1() isa StochasticDiffEqNewtonAdaptiveAlgorithm
    end
end

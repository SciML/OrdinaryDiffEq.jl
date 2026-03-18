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

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence1"
    @time @safetestset "Platen's PL1WM weak second order" begin
        include("weak_convergence/PL1WM.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence2"
    @time @safetestset "Roessler weak SRK Tests" begin
        include("weak_convergence/srk_weak_final.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "W2Ito1WeakConvergence"
    @time @safetestset "Tang & Xiao weak SRK Tests" begin
        include("weak_convergence/W2Ito1.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence3"
    @time @safetestset "Roessler weak SRK (non-diagonal) Tests" begin
        include("weak_convergence/srk_weak_final_non_diagonal.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "SIESMEWeakConvergence"
    @time @safetestset "SIE SME weak Tests" begin
        include("weak_convergence/SIE_SME.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence4"
    @time @safetestset "Weak Stratonovich (non-diagonal) Tests" begin
        include("weak_convergence/weak_strat_non_diagonal.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence5"
    @time @safetestset "Weak Stratonovich Tests" begin
        include("weak_convergence/weak_strat.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakConvergence6"
    @time @safetestset "Roessler weak SRK diagonal Tests" begin
        include("weak_convergence/srk_weak_diagonal_final.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "IRI1WeakConvergence"
    @time @safetestset "IRI1 Weak Convergence Tests" begin
        include("weak_convergence/iri1_weak.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "WeakAdaptiveCPU"
    @time @safetestset "CPU Weak adaptive" begin
        include("adaptive/sde_weak_adaptive.jl")
    end
    @time @safetestset "CPU Weak adaptive step size Brusselator" begin
        include("adaptive/sde_weak_brusselator_adaptive.jl")
    end
end

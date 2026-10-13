using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqHighOrder
        using Test

        @test SRI() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRIW1() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRIW2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRI() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRI2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA1() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA3() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRA() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRA2() isa StochasticDiffEqAdaptiveAlgorithm
    end
    @time @safetestset "SRACache SciMLBase cache hooks" include("sra_cache_hooks_tests.jl")

    @time @safetestset "SRA/SRI resize! stage buffers (issue 4720)" begin
        include("sra_resize_tests.jl")
    end

    @time @safetestset "SRA1 iip k₁ scale allocation" begin
        include("sra1_alloc_tests.jl")
    end

    @time @safetestset "Float32 OOP stages and SRI() Vector OOP" begin
        include("float32_oop_stages_tests.jl")
    end
    @time @safetestset "SRA1 OOP non-diagonal additive noise" begin
        include("sra1_oop_nondiag_tests.jl")
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

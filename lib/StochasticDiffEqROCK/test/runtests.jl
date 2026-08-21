using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqROCK
        using Test

        @test SROCK1() isa StochasticDiffEqAlgorithm
        @test SROCK2() isa StochasticDiffEqAlgorithm
        @test KomBurSROCK2() isa StochasticDiffEqAlgorithm
        @test SROCKC2() isa StochasticDiffEqAlgorithm
        @test SROCKEM() isa StochasticDiffEqAlgorithm
        @test SKSROCK() isa StochasticDiffEqAlgorithm
        @test TangXiaoSROCK2() isa StochasticDiffEqAlgorithm
    end

    @time @safetestset "KomBurSROCK2 non-diagonal noise" begin
        include("kombursrock2_nondiag_tests.jl")
    end

    @time @safetestset "SROCK non-diagonal noise regressions" begin
        include("srock_nondiag_tests.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "SROCKC2WeakConvergence"
    @time @safetestset "SROCKC2 Weak Convergence Tests" begin
        include("weak_convergence/weak_srockc2.jl")
    end
end

if TEST_GROUP == "ALL" || TEST_GROUP == "SROCK2NonDiagonalConvergence"
    @time @safetestset "SROCK2 Non-Diagonal Weak Convergence Tests" begin
        include("weak_convergence/weak_srock2_nondiag.jl")
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

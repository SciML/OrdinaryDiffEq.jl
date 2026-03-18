using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads" begin
        using StochasticDiffEqLowOrder
        using Test

        @test isdefined(StochasticDiffEqLowOrder, :EM)
        @test isdefined(StochasticDiffEqLowOrder, :EulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEM)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :SimplifiedEM)
        @test isdefined(StochasticDiffEqLowOrder, :SplitEM)
        @test isdefined(StochasticDiffEqLowOrder, :RKMil)
        @test isdefined(StochasticDiffEqLowOrder, :RKMilCommute)
        @test isdefined(StochasticDiffEqLowOrder, :PCEuler)
    end
end

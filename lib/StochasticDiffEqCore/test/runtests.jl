using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and types" begin
        using StochasticDiffEqCore
        using Test

        @testset "Module loads" begin
            @test isdefined(StochasticDiffEqCore, :SDEIntegrator)
            @test isdefined(StochasticDiffEqCore, :StochasticDiffEqAlgorithm)
            @test isdefined(StochasticDiffEqCore, :StochasticDiffEqAdaptiveAlgorithm)
            @test isdefined(StochasticDiffEqCore, :StochasticCompositeAlgorithm)
            @test isdefined(StochasticDiffEqCore, :IICommutative)
            @test isdefined(StochasticDiffEqCore, :IILevyArea)
            @test isdefined(StochasticDiffEqCore, :AutoSwitch)
        end

        @testset "Abstract type hierarchy" begin
            @test StochasticDiffEqNewtonAdaptiveAlgorithm <: StochasticDiffEqAdaptiveAlgorithm
            @test StochasticDiffEqNewtonAlgorithm <: StochasticDiffEqAlgorithm
            @test StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm
            @test StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm
        end

        @testset "Iterated integrals" begin
            @test IICommutative() isa IteratedIntegralApprox
            @test IILevyArea() isa IteratedIntegralApprox
        end
    end

    @time @safetestset "_resolve_rng prefers prob.seed over TaskLocalRNG" begin
        # Regression test for the weak-convergence PL1WM failure:
        # When the ensemble layer passes `rng=Random.default_rng()` (TaskLocalRNG)
        # and `prob.seed` is non-zero (typical after `remake(prob, seed=s)` inside
        # a `prob_func`), the integrator must honor `prob.seed` rather than
        # consuming a seed from the TaskLocalRNG. Otherwise two algorithms that
        # differ only in how much randomness they draw per step (e.g. PL1WM vs
        # PL1WMA for additive noise, where PL1WM needs an extra Z process) will
        # see unrelated noise streams from the same ensemble seed assignment.
        using StochasticDiffEqCore
        using Random
        using DiffEqBase
        using Test

        # Fake minimal RODEProblem-like object that just carries a `.seed` field
        # so we can test _resolve_rng in isolation from the full solve pipeline.
        struct FakeRODEProblem <: DiffEqBase.AbstractRODEProblem{Nothing, Nothing, Nothing, Nothing}
            seed::UInt64
        end

        target_seed = UInt64(0x1234567890abcdef)

        prob = FakeRODEProblem(target_seed)
        rng = Random.default_rng()

        # TaskLocalRNG path with prob.seed set → prob.seed wins
        rng_out, _seed, rng_provided = StochasticDiffEqCore._resolve_rng(rng, UInt64(0), prob)
        @test _seed == target_seed
        @test !rng_provided
        @test rng_out isa Random.Xoshiro

        # Explicit seed kwarg takes priority over prob.seed
        rng_out, _seed, rng_provided = StochasticDiffEqCore._resolve_rng(rng, UInt64(42), prob)
        @test _seed == UInt64(42)
        @test !rng_provided

        # prob.seed == 0 → fall back to TaskLocalRNG-derived seed
        prob0 = FakeRODEProblem(UInt64(0))
        rng_out, _seed, rng_provided = StochasticDiffEqCore._resolve_rng(rng, UInt64(0), prob0)
        @test rng_provided
        @test rng_out isa Random.Xoshiro

        # Non-TaskLocalRNG path is unchanged: explicit user RNG is used as-is
        user_rng = Random.Xoshiro(99)
        rng_out, _seed, rng_provided = StochasticDiffEqCore._resolve_rng(user_rng, UInt64(0), prob)
        @test rng_out === user_rng
        @test rng_provided
    end
end

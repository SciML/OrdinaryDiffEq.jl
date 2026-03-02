using StochasticDiffEq, DiffEqNoiseProcess, JumpProcesses, SciMLBase, Random, Test

# Simple SDE: dX = -X dt + 0.1 dW
f(u, p, t) = -u
g(u, p, t) = 0.1
u0 = 1.0
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)

@testset "Integrator RNG Support" begin
    @testset "Default behavior" begin
        sol = solve(prob, EM(); dt = 0.01)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.t) > 1
    end

    @testset "Reproducibility with explicit rng" begin
        rng1 = Xoshiro(42)
        rng2 = Xoshiro(42)
        sol1 = solve(prob, EM(); dt = 0.01, rng = rng1)
        sol2 = solve(prob, EM(); dt = 0.01, rng = rng2)
        @test sol1.u == sol2.u
        @test sol1.t == sol2.t
    end

    @testset "has_rng interface" begin
        integ = init(prob, EM(); dt = 0.01)
        @test SciMLBase.has_rng(integ) == true
    end

    @testset "get_rng returns provided RNG" begin
        rng = Xoshiro(123)
        integ = init(prob, EM(); dt = 0.01, rng)
        @test SciMLBase.get_rng(integ) === rng
    end

    @testset "set_rng! type safety" begin
        integ = init(prob, EM(); dt = 0.01, rng = Xoshiro(42))
        new_rng = Xoshiro(99)
        SciMLBase.set_rng!(integ, new_rng)
        @test SciMLBase.get_rng(integ) === new_rng

        # Mismatched type should error
        @test_throws ArgumentError SciMLBase.set_rng!(integ, MersenneTwister(42))
    end

    @testset "set_rng! syncs framework W.rng" begin
        rng = Xoshiro(42)
        integ = init(prob, EM(); dt = 0.01, rng)
        @test integ.W.rng === rng  # same object initially

        new_rng = Xoshiro(99)
        SciMLBase.set_rng!(integ, new_rng)
        @test integ.W.rng === new_rng  # synced after set_rng!
        @test integ.rng === new_rng
    end

    @testset "set_rng! respects user-provided noise" begin
        user_rng = Xoshiro(77)
        W = WienerProcess(0.0, 0.0; rng = user_rng)
        prob_noise = SDEProblem(f, g, u0, tspan; noise = W)
        integ = init(prob_noise, EM(); dt = 0.01, rng = Xoshiro(42))
        @test integ.user_provided_noise == true

        orig_W_rng = integ.W.rng
        new_rng = Xoshiro(99)
        SciMLBase.set_rng!(integ, new_rng)
        @test integ.rng === new_rng
        @test integ.W.rng === orig_W_rng  # NOT synced
    end

    @testset "reinit! with rng kwarg" begin
        rng = Xoshiro(42)
        integ = init(prob, EM(); dt = 0.01, rng)

        new_rng = Xoshiro(99)
        reinit!(integ, u0; rng = new_rng)
        @test integ.rng === new_rng
        @test integ.W.rng === new_rng
    end

    @testset "Backward compat with seed kwarg" begin
        sol1 = solve(prob, EM(); dt = 0.01, seed = UInt64(42))
        sol2 = solve(prob, EM(); dt = 0.01, seed = UInt64(42))
        @test sol1.u == sol2.u
    end

    @testset "Concrete RNG type parameter" begin
        integ = init(prob, EM(); dt = 0.01, rng = Xoshiro(42))
        # The rng field should be a concrete Xoshiro, not Any
        @test typeof(integ.rng) === Xoshiro
    end

    @testset "Solution seed when rng provided" begin
        # Non-TaskLocalRNG: seed should be UInt64(0)
        sol = solve(prob, EM(); dt = 0.01, rng = Xoshiro(42))
        @test sol.seed == UInt64(0)

        # TaskLocalRNG: seed should be the derived seed (nonzero)
        Random.seed!(12345)
        sol_tlrng = solve(prob, EM(); dt = 0.01, rng = Random.default_rng())
        @test sol_tlrng.seed != UInt64(0)
    end

    @testset "TaskLocalRNG conversion" begin
        integ = init(prob, EM(); dt = 0.01, rng = Random.default_rng())
        @test integ.rng isa Xoshiro
        @test !(integ.rng isa Random.TaskLocalRNG)
    end

    @testset "TaskLocalRNG no duplicate stream" begin
        Random.seed!(12345)
        integ = init(prob, EM(); dt = 0.01, rng = Random.default_rng())
        # The integrator's RNG should produce different values than the task-local RNG.
        # Compare vectors to make a probabilistic flake essentially impossible.
        vals_integ = rand(integ.rng, 10)
        vals_task = rand(10)
        @test vals_integ != vals_task
    end

    @testset "user_provided_noise flag" begin
        integ = init(prob, EM(); dt = 0.01)
        @test integ.user_provided_noise == false

        W = WienerProcess(0.0, 0.0)
        prob_noise = SDEProblem(f, g, u0, tspan; noise = W)
        integ2 = init(prob_noise, EM(); dt = 0.01)
        @test integ2.user_provided_noise == true
    end

    @testset "rng kwarg type validation" begin
        @test_throws ArgumentError solve(prob, EM(); dt = 0.01, rng = 42)
        @test_throws ArgumentError solve(prob, EM(); dt = 0.01, rng = "not an rng")
    end

    @testset "JumpProblem with rng kwarg" begin
        # Setup: zero-drift SDE with a constant-rate jump
        fj!(du, u, p, t) = (du[1] = 0.0)
        gj!(du, u, p, t) = (du[1] = 0.0)
        jump_rate(u, p, t) = 1.0
        jump_affect!(integrator) = (integrator.u[1] += 1.0)
        jump = ConstantRateJump(jump_rate, jump_affect!)

        sde_prob = SDEProblem(fj!, gj!, [0.0], (0.0, 1.0))
        jprob = JumpProblem(sde_prob, Direct(), jump)

        # Solve with explicit rng works
        sol = solve(jprob, EM(); dt = 0.01, rng = Xoshiro(42))
        @test sol.retcode == ReturnCode.Success

        # Repeated solve works (jump state properly reset)
        sol2 = solve(jprob, EM(); dt = 0.01, rng = Xoshiro(99))
        @test sol2.retcode == ReturnCode.Success

        # Note: full JumpProblem reproducibility requires controlling the
        # aggregation RNG too (Phase 4). Here we only verify the rng kwarg
        # does not error or break jump-diffusion solves.
    end

    @testset "VariableRateJump clock refresh with rng kwarg" begin
        # VariableRateJumps use ExtendedJumpArray, whose jump_u (jump clocks)
        # must be refreshed on each init/solve. This tests that
        # reset_jump_problem! is called even when rng is provided.
        fv!(du, u, p, t) = (du[1] = 0.0)
        gv!(du, u, p, t) = (du[1] = 0.0)
        vr_rate(u, p, t) = 0.5
        vr_affect!(integrator) = (integrator.u[1] += 1.0)
        vrjump = VariableRateJump(vr_rate, vr_affect!)

        sde_prob_vr = SDEProblem(fv!, gv!, [0.0], (0.0, 1.0))
        jprob_vr = JumpProblem(sde_prob_vr, Direct(), vrjump)

        integ1 = init(jprob_vr, EM(); dt = 0.01, rng = Xoshiro(42))
        @test integ1.u isa JumpProcesses.ExtendedJumpArray
        jump_u1 = copy(integ1.u.jump_u)

        # Second init on same aliased problem: jump clocks should be refreshed
        integ2 = init(jprob_vr, EM(); dt = 0.01, rng = Xoshiro(99))
        jump_u2 = copy(integ2.u.jump_u)

        @test jump_u1 != jump_u2
    end
end

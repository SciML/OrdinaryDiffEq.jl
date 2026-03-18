using StochasticDiffEq, JumpProcesses, DiffEqBase, Test

# Test that JumpProblem kwargs (tstops, callbacks) are properly forwarded
# to the SDE/RODE solver via merge_problem_kwargs.

# --- Shared setup ---
# Zero-drift, zero-noise SDE so u stays constant (deterministic for testing).
function f!(du, u, p, t)
    return du[1] = 0.0
end

function g!(du, u, p, t)
    return du[1] = 0.0
end

# A jump that never fires (rate = 0).
never_rate(u, p, t) = 0.0
noop_affect!(integrator) = nothing
noop_jump = ConstantRateJump(never_rate, noop_affect!)

@testset "JumpProblem kwarg forwarding" begin

    @testset "tstops passed to solve are forwarded through JumpProblem" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump)
        sol = solve(jprob, EM(); dt = 0.1, tstops = [3.0, 7.0])
        @test 3.0 ∈ sol.t
        @test 7.0 ∈ sol.t
    end

    @testset "callable tstops passed to solve with JumpProblem" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump)
        my_tstops = (p, tspan) -> [2.5, 5.5]
        sol = solve(jprob, EM(); dt = 0.1, tstops = my_tstops)
        @test 2.5 ∈ sol.t
        @test 5.5 ∈ sol.t
    end

    # The JumpProblem constructor natively stores `callback` in jprob.kwargs.
    # These tests verify that merge_problem_kwargs properly forwards it.

    @testset "callback from JumpProblem.kwargs forwarded and fires" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        condition(u, t, integrator) = t == 2.0
        affect!(integrator) = (integrator.u[1] = 42.0)
        cb = DiscreteCallback(condition, affect!)

        jprob = JumpProblem(prob, Direct(), noop_jump; callback = cb)
        sol = solve(jprob, EM(); dt = 0.01, tstops = [2.0])
        @test 2.0 ∈ sol.t
        # Verify callback fired: u was 1.0 before, should be 42.0 after
        idx = findfirst(==(2.0), sol.t)
        @test sol.u[idx + 1][1] == 42.0
    end

    @testset "callback from JumpProblem.kwargs fires exactly once" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        fire_count = Ref(0)
        cb_cond(u, t, integrator) = t == 2.0
        cb_affect!(integrator) = (fire_count[] += 1; integrator.u[1] = 42.0)
        cb = DiscreteCallback(cb_cond, cb_affect!)

        jprob = JumpProblem(prob, Direct(), noop_jump; callback = cb)
        fire_count[] = 0
        sol = solve(jprob, EM(); dt = 0.01, tstops = [2.0])
        @test fire_count[] == 1
    end

    @testset "callback merging between JumpProblem.kwargs and solve kwargs" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        # Callback in jprob.kwargs (stored natively by JumpProblem constructor)
        cb1_cond(u, t, integrator) = t == 2.0
        cb1_affect!(integrator) = (integrator.u[1] = 42.0)
        cb1 = DiscreteCallback(cb1_cond, cb1_affect!)

        jprob = JumpProblem(prob, Direct(), noop_jump; callback = cb1)

        # Callback in solve kwargs
        cb2_cond(u, t, integrator) = t == 4.0
        cb2_affect!(integrator) = (integrator.u[1] = 99.0)
        cb2 = DiscreteCallback(cb2_cond, cb2_affect!)

        sol = solve(jprob, EM(); dt = 0.01, tstops = [2.0, 4.0], callback = cb2)

        # Both callbacks should have fired
        idx2 = findfirst(==(2.0), sol.t)
        @test sol.u[idx2 + 1][1] == 42.0
        idx4 = findfirst(==(4.0), sol.t)
        @test sol.u[idx4 + 1][1] == 99.0
    end

    @testset "merged callbacks each fire exactly once" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        fire_a = Ref(0)
        fire_b = Ref(0)

        cb_a = DiscreteCallback(
            (u, t, integrator) -> t == 2.0,
            integrator -> (fire_a[] += 1)
        )

        jprob = JumpProblem(prob, Direct(), noop_jump; callback = cb_a)

        cb_b = DiscreteCallback(
            (u, t, integrator) -> t == 4.0,
            integrator -> (fire_b[] += 1)
        )

        fire_a[] = 0
        fire_b[] = 0
        sol = solve(jprob, EM(); dt = 0.01, tstops = [2.0, 4.0], callback = cb_b)
        @test fire_a[] == 1
        @test fire_b[] == 1
    end

    @testset "tstops are not duplicated" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump)
        sol = solve(jprob, EM(); dt = 0.1, tstops = [3.0, 7.0])
        @test count(==(3.0), sol.t) == 1
        @test count(==(7.0), sol.t) == 1
    end

    @testset "callback from JumpProblem.kwargs forwarded via init" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        condition(u, t, integrator) = t == 2.0
        affect!(integrator) = (integrator.u[1] = 42.0)
        cb = DiscreteCallback(condition, affect!)

        jprob = JumpProblem(prob, Direct(), noop_jump; callback = cb)
        integrator = init(jprob, EM(); dt = 0.01, tstops = [2.0])
        solve!(integrator)
        idx = findfirst(==(2.0), integrator.sol.t)
        @test integrator.sol.u[idx + 1][1] == 42.0
    end

    # JumpProcesses 9.22+ stores tstops passed to JumpProblem() in jprob.kwargs.
    # These tests verify that tstops from jprob.kwargs are forwarded to the solver.

    @testset "tstops from JumpProblem.kwargs are forwarded via solve" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump; tstops = [3.0, 7.0])
        sol = solve(jprob, EM(); dt = 0.1)
        @test 3.0 ∈ sol.t
        @test 7.0 ∈ sol.t
    end

    @testset "tstops from JumpProblem.kwargs are forwarded via init" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump; tstops = [3.0, 7.0])
        integrator = init(jprob, EM(); dt = 0.1)
        solve!(integrator)
        @test 3.0 ∈ integrator.sol.t
        @test 7.0 ∈ integrator.sol.t
    end

    @testset "tstops from JumpProblem.kwargs are not duplicated" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump; tstops = [3.0, 7.0])
        sol = solve(jprob, EM(); dt = 0.1)
        @test count(==(3.0), sol.t) == 1
        @test count(==(7.0), sol.t) == 1
    end

    @testset "tstops from JumpProblem.kwargs with event that depends on them" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0))

        # Event at t == 2.0 that sets u to 42.0 — only fires if tstops ensures
        # the solver steps exactly at t = 2.0
        condition(u, t, integrator) = t == 2.0
        affect!(integrator) = (integrator.u[1] = 42.0)
        cb = DiscreteCallback(condition, affect!)

        # Both tstops and callback stored in jprob.kwargs — nothing extra to solve
        jprob = JumpProblem(
            prob, Direct(), noop_jump;
            tstops = [2.0], callback = cb
        )
        sol = solve(jprob, EM(); dt = 0.01)
        @test 2.0 ∈ sol.t
        idx = findfirst(==(2.0), sol.t)
        @test sol.u[idx + 1][1] == 42.0
    end

    @testset "callable tstops from JumpProblem.kwargs are forwarded" begin
        prob = SDEProblem(f!, g!, [1.0], (0.0, 5.0), [0.5])

        my_tstops = (p, tspan) -> [p[1] * 4.0]  # tstop at 2.0
        condition(u, t, integrator) = t == 2.0
        affect!(integrator) = (integrator.u[1] = 42.0)
        cb = DiscreteCallback(condition, affect!)

        jprob = JumpProblem(prob, Direct(), noop_jump; tstops = my_tstops, callback = cb)
        sol = solve(jprob, EM(); dt = 0.01)
        @test 2.0 ∈ sol.t
        idx = findfirst(==(2.0), sol.t)
        @test sol.u[idx + 1][1] == 42.0
    end

    # merge_problem_kwargs gives priority to solve kwargs over prob.kwargs for
    # non-callback kwargs (only callbacks get special merging). Verify solve-level
    # tstops take precedence when both sources provide tstops.
    @testset "solve tstops take precedence over JumpProblem.kwargs tstops" begin
        prob = SDEProblem(f!, g!, [10.0], (0.0, 10.0))
        jprob = JumpProblem(prob, Direct(), noop_jump; tstops = [3.0])
        sol = solve(jprob, EM(); dt = 0.1, tstops = [7.0])
        @test 7.0 ∈ sol.t
        @test count(==(7.0), sol.t) == 1
    end

    @testset "tstops forwarding with adaptive SDE solver (SOSRI)" begin
        f_drift!(du, u, p, t) = (du[1] = -0.5 * u[1])
        g_noise!(du, u, p, t) = (du[1] = 0.1 * u[1])
        prob = SDEProblem(f_drift!, g_noise!, [1.0], (0.0, 5.0))
        jprob = JumpProblem(prob, Direct(), noop_jump)
        sol = solve(jprob, SOSRI(); tstops = [1.5, 3.5])
        @test 1.5 ∈ sol.t
        @test 3.5 ∈ sol.t
    end

    @testset "tstops forwarding with RODE solver" begin
        f_rode(u, p, t, W) = 1.01u + 0.87u * W
        prob = RODEProblem(f_rode, 1.0, (0.0, 5.0))
        jprob = JumpProblem(prob, Direct(), noop_jump)
        sol = solve(jprob, RandomEM(); dt = 0.01, tstops = [2.0, 4.0])
        @test 2.0 ∈ sol.t
        @test 4.0 ∈ sol.t
    end
end

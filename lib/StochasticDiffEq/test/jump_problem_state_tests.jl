using StochasticDiffEq, JumpProcesses, SciMLBase, Random, Test

# Under JumpProcesses 10 a JumpProblem stores no RNG: jumps draw from the integrator's
# RNG, and the jump callbacks are re-initialized when an integrator is set up.
# StochasticDiffEq copies a JumpProblem's jump state, or reuses it, according to the
# `alias_jumps` field of the alias specifier.

f!(du, u, p, t) = (du .= 0; nothing)
g!(du, u, p, t) = (du .= 0.1; nothing)
no_noise!(du, u, p, t) = (du .= 0; nothing)
birth_rate(u, p, t) = p[1]
birth!(integrator) = (integrator.u[1] += 1; nothing)

const T_END = 2.0

function birth_problem(; noise = g!, p = [5.0])
    sde = SDEProblem(f!, noise, [0.0], (0.0, T_END), p)
    return JumpProblem(sde, Direct(), ConstantRateJump(birth_rate, birth!))
end

same_path(sol1, sol2) = sol1.t == sol2.t && sol1.u == sol2.u
completed(sol) = SciMLBase.successful_retcode(sol) && sol.t[end] == T_END

# Two solves given the same RNG, or the same seed, run to the end and follow the same path.
function check_replay(solution)
    # Build the options for each solve, so the two solves never share an RNG object.
    for options in (() -> (; rng = Xoshiro(1)), () -> (; seed = 7))
        sol1 = solution(; options()...)
        sol2 = solution(; options()...)
        @test completed(sol1)
        @test completed(sol2)
        @test same_path(sol1, sol2)
    end
    return nothing
end

# The jump aggregator an integrator uses.
function aggregator(integrator)
    return only(
        cb.condition for cb in integrator.opts.callback.discrete_callbacks
            if cb.condition isa JumpProcesses.AbstractSSAJumpAggregator
    )
end

@testset "Jump-diffusion solves are reproducible" begin
    @testset "$(nameof(typeof(alg)))" for (alg, options) in (
            (EM(), (; dt = 0.01)), (SOSRI(), (;)), (ImplicitEM(), (; dt = 0.01)),
        )
        jprob = birth_problem()
        solution(; kwargs...) = solve(jprob, alg; options..., kwargs...)
        check_replay(solution)
        @test !same_path(solution(; seed = 7), solution(; seed = 8))
    end

    @testset "TauLeaping" begin
        rate!(out, u, p, t) = (out[1] = p[1]; nothing)
        change!(du, u, p, t, counts, mark) = (du[1] = counts[1]; nothing)
        jprob = JumpProblem(
            DiscreteProblem([0.0], (0.0, T_END), [5.0]), Direct(),
            RegularJump(rate!, change!, 1)
        )
        solution(; kwargs...) = solve(jprob, TauLeaping(); kwargs...)
        check_replay(solution)
    end
end

@testset "alias_jumps controls copying of jump state" begin
    jprob = birth_problem()
    copy_state = SciMLBase.SDEAliasSpecifier(; alias_jumps = false)
    reuse_state = SciMLBase.SDEAliasSpecifier(; alias_jumps = true)
    @test aggregator(init(jprob, EM(); dt = 0.01, alias = copy_state)) !==
        jprob.discrete_jump_aggregation
    @test aggregator(init(jprob, EM(); dt = 0.01, alias = reuse_state)) ===
        jprob.discrete_jump_aggregation

    # Copied jump state isolates integrators that are alive at the same time.
    expected1 = solve(jprob, EM(); dt = 0.01, rng = Xoshiro(1), alias = copy_state)
    expected2 = solve(jprob, EM(); dt = 0.01, rng = Xoshiro(2), alias = copy_state)
    integrator1 = init(jprob, EM(); dt = 0.01, rng = Xoshiro(1), alias = copy_state)
    integrator2 = init(jprob, EM(); dt = 0.01, rng = Xoshiro(2), alias = copy_state)
    for _ in 1:50
        step!(integrator1)
    end
    solve!(integrator2)
    solve!(integrator1)
    @test completed(expected1)
    @test completed(expected2)
    @test same_path(integrator1.sol, expected1)
    @test same_path(integrator2.sol, expected2)
end

@testset "Mass-action rates follow parameter changes" begin
    birth = MassActionJump([[0 => 1]], [[1 => 1]]; param_idxs = [1])
    sde = SDEProblem(f!, no_noise!, [0.0], (0.0, T_END), [0.0])
    jprob = JumpProblem(sde, Direct(), birth)
    inactive_sol = solve(jprob, EM(); dt = 0.01, seed = 1)
    @test completed(inactive_sol)
    @test inactive_sol.u[end] == [0.0]

    # Remade parameters are used at the next solve.
    active_sol = solve(remake(jprob; p = [50.0]), EM(); dt = 0.01, seed = 1)
    @test completed(active_sol)
    @test active_sol.u[end][1] > 0

    # A parameter change in a callback takes effect after `reset_aggregated_jumps!`.
    state_at_stop = Ref([NaN])
    condition(u, t, integrator) = t == 1.0
    function stop_births!(integrator)
        integrator.p[1] = 0.0
        reset_aggregated_jumps!(integrator)
        state_at_stop[] = copy(integrator.u)
        return nothing
    end
    sol = solve(
        remake(jprob; p = [50.0]), EM(); dt = 0.01, seed = 1, tstops = [1.0],
        callback = DiscreteCallback(condition, stop_births!)
    )
    @test completed(sol)
    @test state_at_stop[][1] > 0
    @test sol.u[end] == state_at_stop[]
end

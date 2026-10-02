using StochasticDiffEq, JumpProcesses, SciMLBase, Test
using Random, StableRNGs
import StochasticDiffEqCore, OrdinaryDiffEqCore

# `init(::JumpProblem, alg)` for StochasticDiffEq algorithms must reach
# StochasticDiffEqCore's JumpProblem-aware initializer, exactly as `solve` does. These
# tests pin the selected methods and check that `init` followed by `solve!` reproduces
# `solve` under the backend's RNG, seed, alias, and callback policies. JumpProcesses'
# OrdinaryDiffEqCore extension defines its own stochastic `__init` method for
# `JumpProblem{IIP, P}`; the Core methods must win without overwriting it.

# ── Fixtures ────────────────────────────────────────────────────────────────

function sde_jump_problem(; seed = 0, variable = false, callback = nothing)
    f!(du, u, p, t) = (du .= -0.1 .* u; nothing)
    g!(du, u, p, t) = (du .= 0.2; nothing)
    prob = SDEProblem(f!, g!, [10.0], (0.0, 0.5), zeros(Int, 2); seed)
    rate(u, p, t) = 10.0
    affect!(integrator) = (integrator.u[1] += 1; nothing)
    jump = variable ? VariableRateJump(rate, affect!) : ConstantRateJump(rate, affect!)
    return JumpProblem(prob, Direct(), jump; callback)
end

function rode_jump_problem(; seed = 0)
    f!(du, u, p, t, W) = (du .= -0.1 .* u .+ 0.2 .* W; nothing)
    prob = RODEProblem(f!, [10.0], (0.0, 0.5), zeros(Int, 2); seed)
    rate(u, p, t) = 10.0
    affect!(integrator) = (integrator.u[1] += 1; nothing)
    return JumpProblem(prob, Direct(), ConstantRateJump(rate, affect!))
end

function regular_jump_problem()
    function regular_rate(out, u, p, t)
        out[1] = (0.1 / 1000.0) * u[1] * u[2]
        out[2] = 0.01u[2]
        return
    end
    dc = [-1.0 0.0; 1.0 -1.0; 0.0 1.0]
    regular_c(du, u, p, t, counts, mark) = (du .= dc * counts; nothing)
    rj = RegularJump(regular_rate, regular_c, 2)
    prob = DiscreteProblem([999.0, 1, 0], (0.0, 25.0))
    return JumpProblem(prob, Direct(), rj)
end

const SDEAlgUnion = Union{
    StochasticDiffEqCore.StochasticDiffEqAlgorithm,
    StochasticDiffEqCore.StochasticDiffEqRODEAlgorithm,
}
const JumpAlgUnion = Union{
    StochasticDiffEqCore.StochasticDiffEqJumpAlgorithm,
    StochasticDiffEqCore.StochasticDiffEqJumpAdaptiveAlgorithm,
}

init_method(prob, alg) = which(SciMLBase.__init, Tuple{typeof(prob), typeof(alg)})

rng_options(kind) =
    kind === :default ? (;) :
    kind === :seed ? (; seed = 456) :
    kind === :zero_seed ? (; seed = 0) :
    kind === :task_rng_seed ? (; rng = Random.default_rng(), seed = 456) :
    kind === :task_rng ? (; rng = Random.default_rng()) :
    (; rng = StableRNG(789), seed = 456)

# Seed metadata the backend records for each option kind, or `nothing` when the seed is
# drawn at random and only route agreement can be checked.
function expected_seed(kind, stored_seed)
    if kind === :explicit_rng_seed
        return UInt64(0)
    elseif kind in (:seed, :task_rng_seed)
        return 456
    elseif stored_seed != 0
        return stored_seed
    else
        return nothing
    end
end

# ── Dispatch ────────────────────────────────────────────────────────────────

@testset "JumpProblem init dispatches to StochasticDiffEqCore" begin
    sde = sde_jump_problem()
    m = init_method(sde, EM())
    @test m.module === StochasticDiffEqCore
    @test m.sig == Tuple{
        typeof(SciMLBase.__init), JumpProblem,
        StochasticDiffEqCore.StochasticDiffEqAlgorithm,
    }

    rode = rode_jump_problem()
    m = init_method(rode, RandomEM())
    @test m.module === StochasticDiffEqCore
    @test m.sig == Tuple{
        typeof(SciMLBase.__init), JumpProblem,
        StochasticDiffEqCore.StochasticDiffEqRODEAlgorithm,
    }

    # The specialized jump-algorithm initializer must still win for RegularJumps.
    regular = regular_jump_problem()
    m = init_method(regular, TauLeaping())
    @test m.module === StochasticDiffEqCore
    @test m.sig == Tuple{typeof(SciMLBase.__init), JumpProblem, JumpAlgUnion}

    # Plain SDE/RODE problems keep the general method.
    m = init_method(sde.prob, EM())
    @test m.module === StochasticDiffEqCore
    @test m.sig == Tuple{
        typeof(SciMLBase.__init), SciMLBase.AbstractRODEProblem, SDEAlgUnion,
    }
end

@testset "No stochastic init ambiguities introduced" begin
    ext = Base.get_extension(JumpProcesses, :JumpProcessesOrdinaryDiffEqCoreExt)
    mods = ext === nothing ?
        (StochasticDiffEqCore, OrdinaryDiffEqCore, JumpProcesses) :
        (StochasticDiffEqCore, OrdinaryDiffEqCore, JumpProcesses, ext)
    relevant(m) = m.name in (:__init, :__solve, :supports_solve_rng)
    ambiguities = filter(
        pair -> any(relevant, pair),
        Test.detect_ambiguities(mods...)
    )
    # JumpProcesses' generic `__init` and OrdinaryDiffEqCore's broad initializer may
    # report a structural overlap independent of this package; no ambiguity may involve
    # a StochasticDiffEqCore method.
    @test all(pair -> all(m -> m.module !== StochasticDiffEqCore, pair), ambiguities)
end

# ── RNG and seed policy: solve versus init → solve! ─────────────────────────

@testset "solve and init use the same RNG policy" begin
    @testset "$kind stored_seed=$stored_seed $options_kind" for
        (kind, alg, make) in (
                (:SDE_constant, EM(), (; seed) -> sde_jump_problem(; seed)),
                (:SDE_variable, EM(), (; seed) -> sde_jump_problem(; seed, variable = true)),
                (:RODE_constant, RandomEM(), (; seed) -> rode_jump_problem(; seed)),
            ),
            stored_seed in (0, 123),
            options_kind in (
                :default, :seed, :zero_seed, :task_rng_seed, :task_rng,
                :explicit_rng_seed,
            )

        # Reseed the task RNG immediately before each entry point: JumpProcesses 9
        # samples from and may reseed the task-local RNG during construction and
        # initialization, and the backend draws random seeds from it.
        Random.seed!(999)
        sol = solve(make(; seed = stored_seed), alg; dt = 0.01, rng_options(options_kind)...)

        Random.seed!(999)
        options = rng_options(options_kind)
        integrator = init(make(; seed = stored_seed), alg; dt = 0.01, options...)
        if options_kind === :explicit_rng_seed
            @test SciMLBase.get_rng(integrator) === options.rng
        else
            @test SciMLBase.get_rng(integrator) isa Xoshiro
        end
        solve!(integrator)

        @test SciMLBase.successful_retcode(sol)
        @test SciMLBase.successful_retcode(integrator.sol)
        @test sol.t == integrator.sol.t
        @test sol.u == integrator.sol.u
        @test sol.seed == integrator.sol.seed
        expected = expected_seed(options_kind, stored_seed)
        if expected !== nothing
            @test sol.seed == expected
        end
    end

    @testset "RegularJump with TauLeaping" begin
        Random.seed!(999)
        sol = solve(regular_jump_problem(), TauLeaping(); dt = 1.0, adaptive = false, seed = 7)
        Random.seed!(999)
        integrator = init(
            regular_jump_problem(), TauLeaping(); dt = 1.0, adaptive = false, seed = 7
        )
        @test SciMLBase.get_rng(integrator) isa Xoshiro
        solve!(integrator)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == integrator.sol.t
        @test sol.u == integrator.sol.u
    end
end

# ── Callback merging ────────────────────────────────────────────────────────

@testset "init preserves callback merging" begin
    @testset "merge_callbacks=$merge_callbacks entry=$entry" for
        merge_callbacks in (true, false), entry in (:solve, :init)

        cb1 = DiscreteCallback(
            (u, t, integrator) -> t == 0.25,
            integrator -> (integrator.p[1] += 1; nothing);
            save_positions = (false, false),
        )
        cb2 = DiscreteCallback(
            (u, t, integrator) -> t == 0.25,
            integrator -> (integrator.p[2] += 1; nothing);
            save_positions = (false, false),
        )
        jprob = sde_jump_problem(; callback = cb1)
        options = (; seed = 123, dt = 0.01, tstops = [0.25], callback = cb2, merge_callbacks)
        sol = if entry === :solve
            solve(jprob, EM(); options...)
        else
            integrator = init(jprob, EM(); options...)
            solve!(integrator)
            integrator.sol
        end
        @test SciMLBase.successful_retcode(sol)
        # The stored callback fires only when merged; the solve-level one always fires.
        @test sol.prob.p == [Int(merge_callbacks), 1]
    end
end

# ── Jump aliasing through the common alias specifier ────────────────────────

# Whether the original problem's aggregator was advanced tells us whether the solver
# aliased it (`alias_jumps = true`) or worked on a copy (`alias_jumps = false`).
function aggregation_advanced(run, jprob)
    before = jprob.discrete_jump_aggregation.next_jump_time
    run()
    return jprob.discrete_jump_aggregation.next_jump_time != before
end

@testset "init honors the backend alias_jumps policy" begin
    @testset "alias_jumps=$alias_jumps entry=$entry" for
        alias_jumps in (true, false, nothing), entry in (:solve, :init)

        alias = alias_jumps === nothing ? nothing :
            SciMLBase.SDEAliasSpecifier(; alias_jumps, alias_u0 = false)
        expected = alias_jumps === nothing ? (Threads.threadid() == 1) : alias_jumps
        jprob = sde_jump_problem()
        if entry === :solve
            advanced = aggregation_advanced(jprob) do
                sol = solve(jprob, EM(); dt = 0.01, seed = 123, alias)
                @test SciMLBase.successful_retcode(sol)
            end
        else
            integrator = nothing
            advanced = aggregation_advanced(jprob) do
                integrator = init(jprob, EM(); dt = 0.01, seed = 123, alias)
                solve!(integrator)
            end
            aliased_callback = any(integrator.opts.callback.discrete_callbacks) do cb
                cb.condition === jprob.discrete_jump_aggregation
            end
            @test aliased_callback == expected
            @test SciMLBase.successful_retcode(integrator.sol)
            # Unrelated alias choices survive: `alias_u0 = false` copies the state.
            alias_jumps === nothing || @test integrator.u !== jprob.prob.u0
        end
        @test advanced == expected
    end
end

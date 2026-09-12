using JumpProcesses, StochasticDiffEqLeaping

struct AdapterRates{M}
    jump::M
end
(r::AdapterRates)(out, u, p, t) = massaction_rates!(out, r.jump, u)
struct AdapterChange{M}
    jump::M
end
(c::AdapterChange)(du, u, p, t, counts, mark) =
    massaction_stoichiometry_mul!(du, c.jump, counts)

function benchmark_massaction(mode, method, n)
    @assert mode in ("native", "adapter")
    @assert method in ("explicit", "implicit")
    maj = MassActionJump(
        fill(0.1, n), [[i => 1] for i in 1:n],
        [[i => -1, mod1(i + 1, n) => 1] for i in 1:n]
    )
    jump = mode == "native" ? maj : RegularJump(AdapterRates(maj), AdapterChange(maj), n)
    prob = JumpProblem(DiscreteProblem(fill(100.0, n), (0.0, 1.0)), PureLeaping(), jump)
    alg = method == "explicit" ? TauLeaping() : ImplicitTauLeaping()
    run() = solve(prob, alg; dt = 0.01, adaptive = false, seed = 123, save_everystep = false)
    first = @timed run()
    @assert successful_retcode(first.value)
    samples = [@timed(run()) for _ in 1:9]
    return println(
        join(
            (
                mode, method, n, first.time, first.compile_time,
                sort([s.time for s in samples])[5], sort([s.bytes for s in samples])[5],
            ), ','
        )
    )
end
benchmark_massaction(ARGS[1], ARGS[2], parse(Int, ARGS[3]))

using JumpProcesses, StochasticDiffEqLeaping, StochasticDiffEqCore, SciMLBase, Test, Random

regular_leaping_algs = (
    TauLeaping(), CaoTauLeaping(),
    ImplicitTauLeaping(), ThetaTrapezoidalTauLeaping(),
)

@testset "Mass-action input for regular leaping solvers" begin
    maj = MassActionJump(
        [0.02, 0.1], [[1 => 2], [2 => 1]],
        [[1 => -2, 2 => 1], [1 => 2, 2 => -1]]
    )
    prob = DiscreteProblem([40.0, 20.0], (0.0, 0.2))
    function rate!(out, u, p, t)
        out[1] = u[1] * max(u[1] - 1, 0) * 0.01
        out[2] = 0.1 * u[2]
        nothing
    end
    function change!(du, u, p, t, counts, mark)
        du[1] = -2 * counts[1] + 2 * counts[2]
        du[2] = counts[1] - counts[2]
        nothing
    end
    rj = RegularJump(rate!, change!, 2)
    for alg in regular_leaping_algs,
            adaptive in (alg isa Union{TauLeaping, CaoTauLeaping} ? (false, true) : (false,))
        @testset "$(nameof(typeof(alg))), adaptive=$adaptive" begin
            actual = JumpProblem(prob, PureLeaping(), maj; rng = Xoshiro(123))
            reference = JumpProblem(prob, PureLeaping(), rj; rng = Xoshiro(123))
            opts = alg isa SimpleTauLeaping ? (; dt = 0.005, seed = 123) :
                (; dt = 0.005, seed = 123, adaptive)
            sol = solve(actual, alg; opts...)
            ref = solve(reference, alg; opts...)
            @test successful_retcode(sol)
            @test actual.regular_jump === nothing
            if !adaptive
                @test sol.t == ref.t
                @test sol.u == ref.u
            end
            @test all(u -> u[1] + 2u[2] ≈ 80, sol.u)
        end
    end
end

@testset "General regular rates remain supported" begin
    rate!(out, u, p, t) = (out[1] = 100 * u[1] / (1 + u[1]))
    function change!(du, u, p, t, counts, mark)
        du[1] = -counts[1]
        du[2] = counts[1]
        nothing
    end
    jp = JumpProblem(
        DiscreteProblem([100.0, 0.0], (0.0, 0.1)),
        PureLeaping(), RegularJump(rate!, change!, 1)
    )
    for alg in regular_leaping_algs
        opts = alg isa SimpleTauLeaping ? (; dt = 0.001, seed = 42) :
            (; dt = 0.001, seed = 42, adaptive = false)
        sol = solve(jp, alg; opts...)
        @test successful_retcode(sol)
        @test sol.u[end][2] > 0
        @test all(u -> sum(u) ≈ 100, sol.u)
    end
end

@testset "Native rates and implicit drift with zero populations" begin
    maj = MassActionJump(
        [2.0, 6.0], [Pair{Int, Int}[], [1 => 3]],
        [[1 => 1], [1 => -3, 2 => 1]]
    )
    for iip in (true, false), alg in regular_leaping_algs
        f = iip ? ((du, u, p, t) -> copyto!(du, u)) : ((u, p, t) -> u)
        prob = DiscreteProblem{iip}(f, [5.0, 0.0], (0.0, 0.02))
        jp = JumpProblem(prob, PureLeaping(), maj)
        data = jump_noise_data(alg, jp, prob.u0, prob.p, 0.0)
        @test data.c === jp.massaction_jump
        @test data.rate_constants == [2.0, 1.0]
        @test all(isfinite, data.rate_constants)
        sol = solve(jp, alg; dt = 0.001, adaptive = false, seed = 123)
        @test successful_retcode(sol)
        @test all(u -> all(isfinite, u), sol.u)
    end
end

@testset "Unsupported mixed jump representations" begin
    maj = MassActionJump([0.1], [[1 => 1]], [[1 => -1]])
    rj = RegularJump(
        (out, u, p, t) -> (out[1] = u[1]),
        (du, u, p, t, counts, mark) -> (du[1] = -counts[1]), 1
    )
    jp = JumpProblem(
        DiscreteProblem([10.0], (0.0, 1.0)), PureLeaping(),
        JumpSet(; massaction_jumps = maj, regular_jumps = rj)
    )
    for alg in regular_leaping_algs
        @test_throws ArgumentError solve(jp, alg; dt = 0.01, adaptive = false)
    end
end

@testset "Out-of-place general rates" begin
    rate(u, p, t) = [100 * u[1] / (1 + u[1])]
    change(u, p, t, counts, mark) = [-counts[1], counts[1]]
    jp = JumpProblem(
        DiscreteProblem([100.0, 0.0], (0.0, 0.1)), PureLeaping(),
        RegularJump(rate, change, 1)
    )
    for alg in regular_leaping_algs
        sol = solve(jp, alg; dt = 0.001, seed = 42, adaptive = false)
        @test successful_retcode(sol)
        @test sol.u[end][2] > 0
        @test all(u -> sum(u) ≈ 100, sol.u)
    end
end

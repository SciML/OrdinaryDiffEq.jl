using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqImplicit
        using SciMLBase: SDEProblem, solve, successful_retcode
        using Test

        @test ImplicitEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitRKMil() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test STrapezoid() isa ImplicitEM
        @test SImplicitMidpoint() isa ImplicitEM
        @test ISSEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ISSEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test SKenCarp() isa StochasticDiffEqNewtonAdaptiveAlgorithm

        f = (du, u, p, t) -> (du .= -u; nothing)
        g = (du, u, p, t) -> (du .= zero(eltype(u)); nothing)
        prob = SDEProblem(f, g, ones(2), (0.0, 1.0))
        sol = solve(prob, ISSEM(); dt = 0.1, adaptive = false)
        @test successful_retcode(sol)
    end

    @time @safetestset "Non-diagonal noise, out-of-place, fixed step (OrdinaryDiffEq.jl#3172)" begin
        using StochasticDiffEqImplicit
        using SciMLBase: SDEProblem, solve, successful_retcode
        using LinearAlgebra: norm
        using Test

        f_oop(u, p, t) = [-u[1] + 0.1 * u[2], -u[2] - 0.1 * u[1]]
        g_oop(u, p, t) = [0.3 0.1; 0.0 0.2]
        f_iip(du, u, p, t) = (du .= f_oop(u, p, t); nothing)
        g_iip(du, u, p, t) = (du .= g_oop(u, p, t); nothing)
        u0 = [1.0, 0.0]
        tspan = (0.0, 1.0)
        prob_oop = SDEProblem(f_oop, g_oop, u0, tspan; noise_rate_prototype = zeros(2, 2))
        prob_iip = SDEProblem(f_iip, g_iip, u0, tspan; noise_rate_prototype = zeros(2, 2))

        # Same seed and fixed step: out-of-place and in-place must agree
        for alg in (ImplicitEM(), ImplicitEulerHeun())
            sol_oop = solve(prob_oop, alg; dt = 0.01, adaptive = false, seed = 7)
            sol_iip = solve(prob_iip, alg; dt = 0.01, adaptive = false, seed = 7)
            @test successful_retcode(sol_oop)
            @test sol_oop.u[end] ≈ sol_iip.u[end] rtol = 1.0e-10
        end

        # With zero noise the scheme reduces to the deterministic trapezoid rule
        g_zero(u, p, t) = zeros(2, 2)
        prob_zero = SDEProblem(f_oop, g_zero, u0, tspan; noise_rate_prototype = zeros(2, 2))
        A = [-1.0 0.1; -0.1 -1.0]
        sol_zero = solve(prob_zero, ImplicitEM(); dt = 0.001, adaptive = false)
        @test norm(sol_zero.u[end] - exp(A) * u0) < 1.0e-3

        # Adaptive stepping for non-diagonal noise in the out-of-place ImplicitEM is not
        # supported yet (the error estimate assumes diagonal noise): see OrdinaryDiffEq.jl#3172
        @test_broken successful_retcode(solve(prob_oop, ImplicitEM(); dt = 0.01, seed = 7))
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqMilstein
        using Test

        @test RKMilGeneral() isa StochasticDiffEqAdaptiveAlgorithm
        @test WangLi3SMil_A() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_B() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_C() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_D() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_E() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_F() isa StochasticDiffEqAlgorithm
    end

    @time @safetestset "RKMilGeneral diagonal vector noise (issue #3863)" begin
        using StochasticDiffEqMilstein
        using StochasticDiffEqMilstein.SciMLBase
        using LinearAlgebra, Test

        # Reproduction from the issue: diagonal-noise vector problem crashed with a
        # DimensionMismatch because _compute_iterated_I returned a full m×m matrix
        # while the diagonal perform_step! branch expects an element-wise m-vector.
        f_lin(u, p, t) = 1.01 .* u
        g_lin(u, p, t) = 0.87 .* u
        prob = SDEProblem(f_lin, g_lin, [0.5, 0.25], (0.0, 1.0))
        sol = solve(prob, RKMilGeneral(), seed = UInt64(7))
        @test SciMLBase.successful_retcode(sol)

        # In-place diagonal path must work too and match the out-of-place result
        # at a fixed dt (identical noise realization).
        f_iip(du, u, p, t) = (du .= 1.01 .* u)
        g_iip(du, u, p, t) = (du .= 0.87 .* u)
        prob_iip = SDEProblem(f_iip, g_iip, [0.5, 0.25], (0.0, 1.0))
        for seed in (UInt64(7), UInt64(42))
            s_oop = solve(prob, RKMilGeneral(); dt = 1 // 2^8, adaptive = false, seed = seed)
            s_iip = solve(prob_iip, RKMilGeneral(); dt = 1 // 2^8, adaptive = false, seed = seed)
            @test SciMLBase.successful_retcode(s_oop)
            @test SciMLBase.successful_retcode(s_iip)
            @test s_oop.u[end] ≈ s_iip.u[end] rtol = 1.0e-10
        end

        # Strong order ≈ 1.0 against the exact geometric-Brownian-motion solution,
        # evaluated on the solver's own realized Brownian path.
        mu = [1.01, -0.5]; sigma = [0.87, 0.3]; u0 = [0.5, 0.25]
        fg(u, p, t) = mu .* u
        gg(u, p, t) = sigma .* u
        gbm = SDEProblem(fg, gg, u0, (0.0, 1.0))
        function strong_error(dt; ntraj = 200)
            total = 0.0
            for k in 1:ntraj
                sol = solve(
                    gbm, RKMilGeneral(); dt = dt, adaptive = false,
                    seed = UInt64(1000 + k), save_noise = true
                )
                WT = sol.W.W[end]
                exact = u0 .* exp.((mu .- sigma .^ 2 ./ 2) .* 1.0 .+ sigma .* WT)
                total += norm(sol.u[end] .- exact)
            end
            return total / ntraj
        end
        e_coarse = strong_error(1 // 2^4)
        e_fine = strong_error(1 // 2^8)
        @test e_fine < e_coarse
        order = log2(e_coarse / e_fine) / log2((1 // 2^4) / (1 // 2^8))
        @test order > 0.8
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

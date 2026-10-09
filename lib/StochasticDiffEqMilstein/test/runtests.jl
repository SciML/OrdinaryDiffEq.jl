using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

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
            s_oop = solve(prob, RKMilGeneral(); dt = 1 // 2^8, adaptive = false, seed)
            s_iip = solve(prob_iip, RKMilGeneral(); dt = 1 // 2^8, adaptive = false, seed)
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
                    gbm, RKMilGeneral(); dt, adaptive = false,
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

    @time @safetestset "RKMilGeneral Ito diagonal update (scalar iip + SA oop)" begin
        using StochasticDiffEqMilstein
        using StochasticDiffEqMilstein.SciMLBase
        using DiffEqNoiseProcess
        using LinearAlgebra, StaticArrays, Test

        seed = UInt64(7)
        A = [
            -1.0 0.1 0.0 0.0
            0.1 -1.0 0.1 0.0
            0.0 0.1 -1.0 0.1
            0.0 0.0 0.1 -1.0
        ]
        B = [0.1 0.05; 0.05 0.1; 0.02 0.03; 0.1 0.0]
        p = (A = A, σ = 0.1, B = B)
        u0 = [1.0, 0.5, 0.25, 0.125]
        f!(du, u, p, t) = mul!(du, p.A, u)
        gmul!(du, u, p, t) = (du .= p.σ .* u)
        gnd!(du, u, p, t) = (du .= p.B .* u)
        f_oop(u, p, t) = p.A * u
        gnd_oop(u, p, t) = p.B .* u

        T = eltype(u0)
        prob_scalar = SDEProblem{true}(
            f!, gmul!, copy(u0), (0.0, 1.0), p;
            noise = WienerProcess(zero(T), zero(T), zero(T))
        )
        sol_scalar = solve(
            prob_scalar, RKMilGeneral(); dt = 1 // 2^4, adaptive = false, seed
        )
        @test SciMLBase.successful_retcode(sol_scalar)
        @test all(isfinite, sol_scalar.u[end])

        # OOP SVector + non-diagonal Ito must not mutate an immutable J.
        u0_sa = SVector{4}(u0)
        B_sa = SMatrix{4, 2}(B)
        p_sa = (A = SMatrix{4, 4}(A), σ = 0.1, B = B_sa)
        prob_sa = SDEProblem{false}(
            f_oop, gnd_oop, u0_sa, (0.0, 1.0), p_sa;
            noise_rate_prototype = zeros(SMatrix{4, 2, Float64})
        )
        sol_sa = solve(
            prob_sa, RKMilGeneral(); dt = 1 // 2^6, adaptive = false, seed
        )
        @test SciMLBase.successful_retcode(sol_sa)
        @test all(isfinite, sol_sa.u[end])
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

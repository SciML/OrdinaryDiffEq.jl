using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads" begin
        using StochasticDiffEqLowOrder
        using Test

        @test isdefined(StochasticDiffEqLowOrder, :EM)
        @test isdefined(StochasticDiffEqLowOrder, :EulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEM)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :SimplifiedEM)
        @test isdefined(StochasticDiffEqLowOrder, :SplitEM)
        @test isdefined(StochasticDiffEqLowOrder, :RKMil)
        @test isdefined(StochasticDiffEqLowOrder, :RKMilCommute)
        @test isdefined(StochasticDiffEqLowOrder, :PCEuler)
    end

    @time @safetestset "RKMil OOP adaptive error estimate (ggprime)" begin
        using StochasticDiffEqLowOrder
        using StochasticDiffEqLowOrder.SciMLBase
        using DiffEqNoiseProcess
        using LinearAlgebra, Test

        # OOP RKMilConstantCache previously used the drift f(K) in the noise
        # error estimate instead of ggprime = (g(utilde)-L)/sqdt (as IIP does),
        # making En O(1/sqrt(dt)) too large and forcing ~1000x more steps.
        seed = UInt64(7)
        A = [
            -1.0 0.1 0.0 0.0
            0.1 -1.0 0.1 0.0
            0.0 0.1 -1.0 0.1
            0.0 0.0 0.1 -1.0
        ]
        p = (A = A, σ = 0.1)
        u0 = [1.0, 0.5, 0.25, 0.125]
        f_oop(u, p, t) = p.A * u
        g_oop(u, p, t) = p.σ .* u
        f_iip!(du, u, p, t) = mul!(du, p.A, u)
        g_iip!(du, u, p, t) = (du .= p.σ .* u)

        for noise in (:diag, :scalar)
            kw = noise === :scalar ? (; noise = WienerProcess(0.0, 0.0, 0.0)) : (;)
            prob_oop = SDEProblem{false}(f_oop, g_oop, copy(u0), (0.0, 1.0), p; kw...)
            prob_iip = SDEProblem{true}(f_iip!, g_iip!, copy(u0), (0.0, 1.0), p; kw...)
            sol_oop = solve(prob_oop, RKMil(); seed)
            sol_iip = solve(prob_iip, RKMil(); seed)
            @test SciMLBase.successful_retcode(sol_oop)
            @test SciMLBase.successful_retcode(sol_iip)
            @test sol_oop.stats.naccept <= 4 * sol_iip.stats.naccept
            @test maximum(abs.(sol_oop.u[end] .- sol_iip.u[end])) < 1.0e-3
        end

        # Fixed-dt strong order ≈ 1: error should fall roughly ×4 when dt → dt/4.
        mu = 1.01
        sigma = 0.87
        f_lin(u, p, t) = mu * u
        g_lin(u, p, t) = sigma * u
        analytic(u0, p, t, W) = u0 * exp((mu - sigma^2 / 2) * t + sigma * W)
        gbm = SDEProblem(
            SDEFunction(f_lin, g_lin; analytic = analytic), 0.5, (0.0, 1.0)
        )
        function mean_strong_error(dt; ntraj = 200)
            total = 0.0
            for k in 1:ntraj
                sol = solve(
                    gbm, RKMil(); dt, adaptive = false, seed = UInt64(1000 + k)
                )
                total += abs(sol.u[end] - sol.u_analytic[end])
            end
            return total / ntraj
        end
        dt_coarse = 1 // 2^4
        dt_fine = 1 // 2^8
        e_coarse = mean_strong_error(dt_coarse)
        e_mid = mean_strong_error(1 // 2^6)
        e_fine = mean_strong_error(dt_fine)
        @test e_mid < e_coarse
        @test e_fine < e_mid
        order = log(e_coarse / e_fine) / log(Float64(dt_coarse / dt_fine))
        @test order > 0.75
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

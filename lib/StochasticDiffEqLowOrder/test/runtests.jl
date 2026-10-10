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

    @time @safetestset "RKMilCommute Ito diagonal update (scalar iip + nondiag alloc)" begin
        using StochasticDiffEqLowOrder
        using StochasticDiffEqLowOrder.SciMLBase
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
        gmul_oop(u, p, t) = p.σ .* u

        T = eltype(u0)
        prob_scalar = SDEProblem{true}(
            f!, gmul!, copy(u0), (0.0, 1.0), p;
            noise = WienerProcess(zero(T), zero(T), zero(T))
        )
        sol_scalar = solve(
            prob_scalar, RKMilCommute(); dt = 1 // 2^4, adaptive = false, seed
        )
        @test SciMLBase.successful_retcode(sol_scalar)
        @test all(isfinite, sol_scalar.u[end])

        # Scalar IIP must match scalar OOP bitwise at fixed dt.
        prob_scalar_oop = SDEProblem{false}(
            f_oop, gmul_oop, copy(u0), (0.0, 1.0), p;
            noise = WienerProcess(zero(T), zero(T), zero(T))
        )
        sol_scalar_oop = solve(
            prob_scalar_oop, RKMilCommute(); dt = 1 // 2^4, adaptive = false, seed
        )
        @test SciMLBase.successful_retcode(sol_scalar_oop)
        @test sol_scalar.u[end] == sol_scalar_oop.u[end]

        prob_nd = SDEProblem{true}(
            f!, gnd!, copy(u0), (0.0, 1.0), p;
            noise_rate_prototype = zeros(eltype(u0), length(u0), 2)
        )
        integ = init(
            prob_nd, RKMilCommute(); dt = 1 // 2^6, adaptive = false,
            save_everystep = false, seed
        )
        for _ in 1:5
            step!(integ) # warmup compile / noise setup
        end
        allocs = @allocated step!(integ)
        @test allocs < 48

        # OOP SVector + non-diagonal Ito: J is an SMatrix and must not be mutated.
        u0_sa = SVector{4}(u0)
        B_sa = SMatrix{4, 2}(B)
        p_sa = (A = SMatrix{4, 4}(A), σ = 0.1, B = B_sa)
        prob_sa = SDEProblem{false}(
            f_oop, gnd_oop, u0_sa, (0.0, 1.0), p_sa;
            noise_rate_prototype = zeros(SMatrix{4, 2, Float64})
        )
        sol_sa_fixed = solve(
            prob_sa, RKMilCommute(); dt = 1 // 2^6, adaptive = false, seed
        )
        sol_sa_adapt = solve(prob_sa, RKMilCommute(); seed)
        @test SciMLBase.successful_retcode(sol_sa_fixed)
        @test SciMLBase.successful_retcode(sol_sa_adapt)
        @test all(isfinite, sol_sa_fixed.u[end])
        @test all(isfinite, sol_sa_adapt.u[end])
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end

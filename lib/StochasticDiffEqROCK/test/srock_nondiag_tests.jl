using StochasticDiffEqROCK, DiffEqNoiseProcess, Test, Random
import SciMLBase

# Non-square noise (n = 2 states, m = 1 Brownian motion) used to crash with
# MethodError / DimensionMismatch because a length-one W.dW shortcut routed these
# problems into the diagonal-noise code path, see SciML/OrdinaryDiffEq.jl#3170.
@testset "m=1 non-diagonal noise does not crash" begin
    n, m = 2, 1
    f_oop(u, p, t) = zeros(n)
    g_oop(u, p, t) = ones(n, m)
    f_iip(du, u, p, t) = (du .= 0.0)
    g_iip(du, u, p, t) = (du .= 1.0)
    u0 = zeros(n)
    nrp = ones(n, m)
    prob_oop = SDEProblem(f_oop, g_oop, u0, (0.0, 1.0); noise_rate_prototype = nrp)
    prob_iip = SDEProblem(f_iip, g_iip, u0, (0.0, 1.0); noise_rate_prototype = nrp)

    @testset "$(nameof(typeof(alg)))" for alg in (
            SROCK2(), SROCKEM(), SROCKC2(), KomBurSROCK2(), TangXiaoSROCK2(),
        )
        for prob in (prob_oop, prob_iip)
            Random.seed!(100)
            sol = solve(prob, alg, dt = 1 / 2^4)
            @test SciMLBase.successful_retcode(sol)
            @test all(isfinite, sol.u[end])
        end
    end
end

# The Rademacher variables of the general-noise paths are drawn from the RNG of
# the noise process. A NoiseWrapper keeps that RNG on its wrapped source, which
# may itself be a wrapper (see SciML/OrdinaryDiffEq.jl#3188), and a NoiseGrid
# replays a recorded path and carries no RNG at all.
@testset "supplied noise processes (non-diagonal noise)" begin
    f_iip(du, u, p, t) = (du .= -0.5 .* u)
    function g_iip(du, u, p, t)
        du[1, 1] = 0.1 * u[1]
        du[1, 2] = 0.05 * u[2]
        du[2, 1] = 0.05 * u[1]
        du[2, 2] = 0.1 * u[2]
        return nothing
    end
    u0 = [1.0, 1.0]
    tspan = (0.0, 1.0)
    dt = 1 / 2^6

    @testset "$(nameof(typeof(alg)))" for alg in (SROCK2(), KomBurSROCK2())
        Random.seed!(100)
        prob = SDEProblem(f_iip, g_iip, u0, tspan; noise_rate_prototype = zeros(2, 2))
        sol = solve(prob, alg; dt, adaptive = false, save_noise = true)
        @test SciMLBase.successful_retcode(sol)

        for noise in (NoiseWrapper(sol.W), NoiseWrapper(NoiseWrapper(sol.W)))
            prob_wrapped = SDEProblem(
                f_iip, g_iip, u0, tspan;
                noise, noise_rate_prototype = zeros(2, 2)
            )
            sol_wrapped = solve(prob_wrapped, alg; dt, adaptive = false)
            @test SciMLBase.successful_retcode(sol_wrapped)
            @test all(isfinite, sol_wrapped.u[end])
        end

        ts = collect(tspan[1]:dt:tspan[2])
        Ws = Vector{Vector{Float64}}(undef, length(ts))
        Ws[1] = zeros(2)
        for i in 2:length(ts)
            Ws[i] = Ws[i - 1] .+ sqrt(dt) .* randn(2)
        end
        prob_grid = SDEProblem(
            f_iip, g_iip, u0, tspan;
            noise = NoiseGrid(ts, Ws), noise_rate_prototype = zeros(2, 2)
        )
        sol_grid = solve(prob_grid, alg; dt, adaptive = false)
        @test SciMLBase.successful_retcode(sol_grid)
        @test all(isfinite, sol_grid.u[end])
    end
end

@testset "definite assignment paths" begin
    f_oop(u, p, t) = -u
    g_diag_oop(u, p, t) = 0.1u
    g_full_oop(u, p, t) = [0.1u[1] 0.05u[2]; 0.05u[1] 0.1u[2]]
    f_iip(du, u, p, t) = (du .= -u)
    g_diag_iip(du, u, p, t) = (du .= 0.1 .* u)
    function g_full_iip(du, u, p, t)
        du[1, 1] = 0.1u[1]
        du[1, 2] = 0.05u[2]
        du[2, 1] = 0.05u[1]
        du[2, 2] = 0.1u[2]
        return nothing
    end

    u0 = [1.0, 0.5]
    tspan = (0.0, 0.125)
    nrp = zeros(2, 2)
    problems = (
        SDEProblem(f_oop, g_diag_oop, u0, tspan),
        SDEProblem(f_iip, g_diag_iip, u0, tspan),
        SDEProblem(f_oop, g_full_oop, u0, tspan; noise_rate_prototype = nrp),
        SDEProblem(f_iip, g_full_iip, u0, tspan; noise_rate_prototype = nrp),
    )

    for alg in (SROCK1(), SKSROCK(post_processing = true), TangXiaoSROCK2())
        @testset "$(nameof(typeof(alg)))" begin
            for (i, prob) in enumerate(problems)
                Random.seed!(100 + i)
                sol = solve(prob, alg; dt = 1 / 64, adaptive = false)
                @test SciMLBase.successful_retcode(sol)
                @test all(isfinite, sol.u[end])
            end
        end
    end
end

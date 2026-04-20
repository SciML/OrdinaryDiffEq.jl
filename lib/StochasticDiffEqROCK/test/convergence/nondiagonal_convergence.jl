using StochasticDiffEqROCK, DiffEqDevTools, DiffEqNoiseProcess, Test, Random

# Non-diagonal SDE: 2-component system driven by 1 Brownian motion (n=2, m=1)
# du_i = μ u_i dt + σ u_i dW    (same scalar Brownian motion for both components)
# Weak solution: E[u_i(T)] = u_i(0) exp(μ T)

const μ_test = -0.5
const σ_test = 0.1

f_nd_oop(u, p, t) = μ_test .* u
g_nd_oop(u, p, t) = σ_test .* reshape(u, length(u), 1)

f_nd_iip!(du, u, p, t) = (du .= μ_test .* u)
function g_nd_iip!(du, u, p, t)
    for i in axes(du, 1)
        du[i, 1] = σ_test * u[i]
    end
end

# Analytic (Itô): u_i(t) = u_i(0) exp((μ - σ²/2)t + σ W(t))
analytic_nd(u0, p, t, W) = u0 .* exp.((μ_test - σ_test^2 / 2) * t .+ σ_test .* W[1])

u0 = [1.0, 1.0]
tspan = (0.0, 1.0)

prob_oop = SDEProblem(
    SDEFunction(f_nd_oop, g_nd_oop; analytic = analytic_nd),
    u0, tspan;
    noise_rate_prototype = zeros(2, 1)
)
prob_iip = SDEProblem(
    SDEFunction(f_nd_iip!, g_nd_iip!; analytic = analytic_nd),
    u0, tspan;
    noise_rate_prototype = zeros(2, 1)
)

Random.seed!(100)
dts = 1 .// 2 .^ (6:-1:2)

@testset "SROCK2 non-diagonal OOP weak convergence (order ~2)" begin
    sim = test_convergence(dts, prob_oop, SROCK2(), trajectories = Int(5e4),
                           save_everystep = false, weak_timeseries_errors = false)
    @test abs(sim.𝒪est[:weak_final] - 2.0) < 0.5
end

@testset "SROCK2 non-diagonal IIP weak convergence (order ~2)" begin
    sim = test_convergence(dts, prob_iip, SROCK2(), trajectories = Int(5e4),
                           save_everystep = false, weak_timeseries_errors = false)
    @test abs(sim.𝒪est[:weak_final] - 2.0) < 0.5
end

# KomBurSROCK2 smoke test: non-diagonal noise must not crash with DimensionMismatch
@testset "KomBurSROCK2 non-diagonal IIP does not crash" begin
    prob_smoke = SDEProblem(f_nd_iip!, g_nd_iip!, u0, tspan;
                            noise_rate_prototype = zeros(2, 1))
    sol = solve(prob_smoke, KomBurSROCK2(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
end

# Smoke test: NoiseWrapper must not crash
@testset "SROCK2 NoiseWrapper does not crash" begin
    prob_base = SDEProblem(f_nd_iip!, g_nd_iip!, u0, tspan;
                           noise_rate_prototype = zeros(2, 1))
    sol_base = solve(prob_base, SROCK2(), dt = 0.01, save_noise = true)
    W_wrap = NoiseWrapper(sol_base.W)
    prob_wrap = SDEProblem(f_nd_iip!, g_nd_iip!, u0, tspan;
                           noise = W_wrap, noise_rate_prototype = zeros(2, 1))
    sol_wrap = solve(prob_wrap, SROCK2(), dt = 0.01)
    @test sol_wrap.retcode == ReturnCode.Success
end

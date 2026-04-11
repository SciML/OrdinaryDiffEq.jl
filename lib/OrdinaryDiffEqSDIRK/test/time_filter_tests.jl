using OrdinaryDiffEqSDIRK, ODEProblemLibrary, DiffEqDevTools
using SciMLBase: ReturnCode
using Test, Random
Random.seed!(100)

total_time = @elapsed begin

testTol = 0.2

# ──────────────────────────────────────────────────────────────────
# 1. Convergence order — proves mathematical correctness
# ──────────────────────────────────────────────────────────────────

@testset "Convergence order ($(["OOP", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]
    dts = 1 .// 2 .^ (9:-1:5)

    sim_filtered = test_convergence(dts, prob, ImplicitEuler(time_filter = true))
    @test sim_filtered.𝒪est[:final] ≈ 2 atol = testTol

    sim_unfiltered = test_convergence(dts, prob, ImplicitEuler())
    @test sim_unfiltered.𝒪est[:final] ≈ 1 atol = testTol
end

# ──────────────────────────────────────────────────────────────────
# 2. Quantitative accuracy on dy/dt = -y
# ──────────────────────────────────────────────────────────────────

@testset "Analytical accuracy dy/dt = -y" begin
    f(u, p, t) = -u
    prob = ODEProblem(f, 1.0, (0.0, 1.0))
    exact = exp(-1.0)

    for dt in [1 // 64, 1 // 128]
        sol_plain = solve(prob, ImplicitEuler(); dt, adaptive = false)
        sol_filt = solve(prob, ImplicitEuler(time_filter = true); dt, adaptive = false)

        err_plain = abs(sol_plain[end] - exact)
        err_filt = abs(sol_filt[end] - exact)

        # Filtered error must be significantly smaller
        @test err_filt < err_plain
    end

    # Check O(dt²) scaling: halving dt should reduce error by ~4×
    dt_coarse = 1 // 32
    dt_fine = 1 // 64
    sol_c = solve(prob, ImplicitEuler(time_filter = true); dt = dt_coarse, adaptive = false)
    sol_f = solve(prob, ImplicitEuler(time_filter = true); dt = dt_fine, adaptive = false)
    err_c = abs(sol_c[end] - exact)
    err_f = abs(sol_f[end] - exact)
    ratio = err_c / err_f
    # For 2nd-order, ratio ≈ 4; allow generous tolerance
    @test 2.5 < ratio < 6.0
end

# ──────────────────────────────────────────────────────────────────
# 3. Stiff system — stability check
# ──────────────────────────────────────────────────────────────────

@testset "Stiff system stability" begin
    # Stiff linear system: λ = -1000
    f(u, p, t) = -1000.0 * u
    prob = ODEProblem(f, 1.0, (0.0, 1.0))

    sol = solve(prob, ImplicitEuler(time_filter = true); dt = 1 // 100, adaptive = false)
    @test isfinite(sol[end])
    @test abs(sol[end]) < 1e-3  # solution decays to ~0

    # In-place version
    f_ip!(du, u, p, t) = (du .= -1000.0 .* u)
    prob_ip = ODEProblem(f_ip!, [1.0], (0.0, 1.0))
    sol_ip = solve(
        prob_ip, ImplicitEuler(time_filter = true); dt = 1 // 100, adaptive = false)
    @test all(isfinite, sol_ip[end])
    @test all(x -> abs(x) < 1e-3, sol_ip[end])
end

# ──────────────────────────────────────────────────────────────────
# 4. Adaptive timestepping with filter
# ──────────────────────────────────────────────────────────────────

@testset "Adaptive timestepping" begin
    prob = ODEProblemLibrary.prob_ode_linear

    sol_filt = solve(prob, ImplicitEuler(time_filter = true))
    @test sol_filt.retcode == ReturnCode.Success
    @test isfinite(sol_filt[end])

    # In-place
    prob_ip = ODEProblemLibrary.prob_ode_2Dlinear
    sol_ip = solve(prob_ip, ImplicitEuler(time_filter = true))
    @test sol_ip.retcode == ReturnCode.Success
    @test all(isfinite, sol_ip[end])
end

# ──────────────────────────────────────────────────────────────────
# 5. Edge cases
# ──────────────────────────────────────────────────────────────────

@testset "Edge cases" begin
    # Single step: filter cannot activate
    f(u, p, t) = -u
    prob1 = ODEProblem(f, 1.0, (0.0, 0.1))
    sol_f1 = solve(prob1, ImplicitEuler(time_filter = true); dt = 0.1, adaptive = false)
    sol_p1 = solve(prob1, ImplicitEuler(); dt = 0.1, adaptive = false)
    @test sol_f1[end] ≈ sol_p1[end] atol = 1e-14

    # Zero RHS: solution must stay at u0
    f_zero(u, p, t) = zero(u)
    prob_z = ODEProblem(f_zero, 42.0, (0.0, 1.0))
    sol_z = solve(prob_z, ImplicitEuler(time_filter = true); dt = 1 // 10, adaptive = false)
    @test all(u -> u ≈ 42.0, sol_z.u)

    # Very large dt: should not crash
    prob_big = ODEProblem(f, 1.0, (0.0, 1.0))
    sol_big = solve(
        prob_big, ImplicitEuler(time_filter = true); dt = 0.5, adaptive = false)
    @test isfinite(sol_big[end])
end

end  # @elapsed

println("\n✅ All time-filter tests passed.")
@info "Total wall-clock time" total_time_seconds = round(total_time; digits = 2)

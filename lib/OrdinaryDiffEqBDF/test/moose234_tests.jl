# MOOSE234 tests — Variable-stepsize, variable-order (2/3/4) embedded method
# DeCaria et al., arXiv:1810.06670v1
#
# MOOSE234 is a VSVO method (like FBDF). The order varies adaptively, so
# test_convergence with fixed dt is not applicable. Tests verify:
#   1. Solver completes without error (OOP + IIP)
#   2. Solution accuracy against analytic solutions
#   3. Adaptive order switching (startup order progression)
#   4. Stiff problem stability (Van der Pol)
#   5. Static array (SVector) and scalar support
#   6. reinit! correctness

using OrdinaryDiffEqBDF, ODEProblemLibrary, DiffEqDevTools, Test
using StaticArrays

# ────────────────────────────────────────────────────────────────────────────
# 1. Basic convergence: solver completes without error
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 basic solve ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    @test_nowarn solve(prob, MOOSE234())

    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
end

# ────────────────────────────────────────────────────────────────────────────
# 2. Accuracy: scalar exponential decay  dy/dt = -y, y(0) = 1
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 accuracy (scalar exponential)" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    exact = exp(-1.0)

    sol_mod = solve(prob, MOOSE234(), abstol = 1e-6, reltol = 1e-6)
    @test isapprox(sol_mod[end], exact, rtol = 1e-4)

    sol_tight = solve(prob, MOOSE234(), abstol = 1e-10, reltol = 1e-10)
    @test isapprox(sol_tight[end], exact, rtol = 1e-6)

    # Tighter tolerance should yield strictly smaller error
    @test abs(sol_tight[end] - exact) < abs(sol_mod[end] - exact)
end

# ────────────────────────────────────────────────────────────────────────────
# 3. Accuracy: in-place 2D linear system  du/dt = p*u
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 accuracy (in-place 2D)" begin
    fiip = (du, u, p, t) -> du .= p .* u
    u0 = [1.0, 2.0]
    prob = ODEProblem(fiip, u0, (0.0, 1.0), -1.0)

    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
    @test isapprox(sol[end], exp(-1.0) .* u0, rtol = 1e-2)
end

# ────────────────────────────────────────────────────────────────────────────
# 4. MOOSE234 vs FBDF accuracy comparison
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 vs FBDF accuracy" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    exact = exp(-1.0)

    sol_moose = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    sol_fbdf = solve(prob, FBDF(), abstol = 1e-8, reltol = 1e-8)

    err_moose = abs(sol_moose[end] - exact)
    err_fbdf = abs(sol_fbdf[end] - exact)

    # Both should achieve sub-tolerance accuracy
    @test err_moose < 1e-6
    @test err_fbdf < 1e-6
end

# ────────────────────────────────────────────────────────────────────────────
# 5. Startup order progression
#    - Step 1: only BDF1/BDF2 history → order 2
#    - Step 2+: BDF3 becomes available → order ≤ 3
#    - Step 3+: FBDF4 filter available → order ≤ 4
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 startup order progression" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 5.0))
    integ = init(prob, MOOSE234(), dt = 0.01)

    # Before any step: order 2, no history
    @test integ.cache.order == 2
    @test integ.cache.iters_from_event == 0

    # Step 1: limited history
    step!(integ)
    @test integ.cache.iters_from_event >= 1
    # Must still be within valid order range
    @test 2 <= integ.cache.order <= 4

    # Step 2
    step!(integ)
    @test integ.cache.iters_from_event >= 2

    # After enough steps, full order range (2–4) should be accessible
    for _ in 1:10
        step!(integ)
    end
    @test integ.cache.iters_from_event >= 4
end

# ────────────────────────────────────────────────────────────────────────────
# 6. Adaptive VSVO: verify the solver actually switches orders
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 adaptive order switching" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0))
    integ = init(prob, MOOSE234(), dt = 0.01)

    orders_seen = Set{Int}()
    for _ in 1:200
        step!(integ)
        push!(orders_seen, integ.cache.order)
        integ.sol.retcode == ReturnCode.Failure && break
    end

    # Should have used at least two different orders during integration
    @test length(orders_seen) >= 2
    # All orders must be in the valid range
    @test all(o -> 2 <= o <= 4, orders_seen)
end

# ────────────────────────────────────────────────────────────────────────────
# 7. Van der Pol oscillator (stiff system)
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 Van der Pol (stiff)" begin
    function vdp!(du, u, p, t)
        μ = p[1]
        du[1] = u[2]
        du[2] = μ * (1 - u[1]^2) * u[2] - u[1]
    end

    # Moderate stiffness (μ = 100)
    prob_mod = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [100.0])
    sol_mod = solve(prob_mod, MOOSE234(), abstol = 1e-6, reltol = 1e-6)
    @test sol_mod.retcode == ReturnCode.Success
    @test sol_mod.t[end] == 6.3

    # High stiffness (μ = 1000) — the paper's benchmark
    prob_stiff = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [1000.0])
    sol_stiff = solve(prob_stiff, MOOSE234(), abstol = 1e-4, reltol = 1e-4)
    @test sol_stiff.retcode == ReturnCode.Success
    @test sol_stiff.t[end] == 6.3
end

# ────────────────────────────────────────────────────────────────────────────
# 8. ROBER stiff chemical kinetics (3-species DAE-like ODE)
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 ROBER stiff system" begin
    function rober!(du, u, p, t)
        y1, y2, y3 = u
        du[1] = -0.04 * y1 + 1e4 * y2 * y3
        du[2] = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2
        du[3] = 3e7 * y2^2
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1e5))
    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == 1e5
    # Conservation: y1 + y2 + y3 = 1
    @test isapprox(sum(sol[end]), 1.0, atol = 1e-6)
end

# ────────────────────────────────────────────────────────────────────────────
# 9. Static Array (SVector) support
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 Static Array (SVector)" begin
    f_oop(u, p, t) = -0.5 * u
    u0_sv = SVector(1.0, 2.0)
    prob_sv = ODEProblem(f_oop, u0_sv, (0.0, 1.0))

    sol = solve(prob_sv, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.u[end] isa SVector
    @test isapprox(sol.u[end], exp(-0.5) * u0_sv, rtol = 1e-3)
end

# ────────────────────────────────────────────────────────────────────────────
# 10. Scalar ODE
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 scalar" begin
    prob_scalar = ODEProblem((u, p, t) -> -0.5 * u, 1.0, (0.0, 1.0))

    sol = solve(prob_scalar, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.u[end] isa Number
    @test isapprox(sol.u[end], exp(-0.5), rtol = 1e-5)
end

# ────────────────────────────────────────────────────────────────────────────
# 11. reinit! correctness
# ────────────────────────────────────────────────────────────────────────────
@testset "MOOSE234 reinit" begin
    foop = (u, p, t) -> p * u
    proboop = ODEProblem(foop, ones(2), (0.0, 10.0), -1.0)

    fiip = (du, u, p, t) -> du .= p .* u
    probiip = ODEProblem(fiip, ones(2), (0.0, 10.0), -1.0)

    for prob in [proboop, probiip]
        integ = init(prob, MOOSE234())
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.sol.t[end] == 10.0

        reinit!(integ, prob.u0)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.sol.t[end] == 10.0
    end
end

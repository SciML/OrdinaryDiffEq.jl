using OrdinaryDiffEqRosenbrock, Test, LinearAlgebra
using OrdinaryDiffEqRosenbrock: JacReuseState
using SciMLBase

# ============================================================================
# Test problems
# ============================================================================

# Van der Pol oscillator (stiff, mu=1e3)
function vanderpol!(du, u, p, t)
    mu = p[1]
    du[1] = u[2]
    return du[2] = mu * (1 - u[1]^2) * u[2] - u[1]
end
vdp_prob = ODEProblem(vanderpol!, [2.0, 0.0], (0.0, 0.5), [1.0e3])

# ROBER chemical kinetics (stiff)
function rober!(du, u, p, t)
    y1, y2, y3 = u
    du[1] = -0.04 * y1 + 1.0e4 * y2 * y3
    du[2] = 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2
    return du[3] = 3.0e7 * y2^2
end
rober_prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e2))

# Simple exponential decay for convergence tests
function expdecay!(du, u, p, t)
    return du .= -u
end
expdecay_prob = ODEProblem(expdecay!, [1.0, 2.0, 3.0], (0.0, 1.0))

# 2D linear ODE for convergence order tests: du/dt = A*u
# u(t) = exp(A*t)*u0, exact solution computable via matrix exponential
function linear2d!(du, u, p, t)
    du[1] = -2.0 * u[1] + u[2]
    return du[2] = u[1] - 2.0 * u[2]
end
u0_2d = [1.0, 0.0]
tspan_2d = (0.0, 1.0)
A_2d = [-2.0 1.0; 1.0 -2.0]
exact_2d = exp(A_2d * tspan_2d[2]) * u0_2d
linear2d_prob = ODEProblem(linear2d!, u0_2d, tspan_2d)

"""
    estimate_order(prob, alg, dts, exact)

Estimate convergence order by solving at multiple step sizes and
fitting log(error) vs log(dt) via least-squares regression.
"""
function estimate_order(prob, alg, dts, exact)
    errors = Float64[]
    for dt in dts
        sol = solve(prob, alg, dt = dt, adaptive = false)
        push!(errors, norm(sol.u[end] - exact))
    end
    # Linear regression: log(error) = order * log(dt) + c
    log_dts = log.(dts)
    log_errs = log.(errors)
    n = length(dts)
    slope = (n * sum(log_dts .* log_errs) - sum(log_dts) * sum(log_errs)) /
        (n * sum(log_dts .^ 2) - sum(log_dts)^2)
    return slope
end

# Mass matrix DAE: M*u' = f(u)
# x' = -x + y, 0 = x - y  =>  y = x, x' = 0  =>  x(t) = y(t) = 1
function mass_dae!(du, u, p, t)
    du[1] = -u[1] + u[2]
    return du[2] = u[1] - u[2]
end
M = [1.0 0.0; 0.0 0.0]
mass_f = ODEFunction(mass_dae!, mass_matrix = M)
mass_prob = ODEProblem(mass_f, [1.0, 1.0], (0.0, 1.0))

# ============================================================================
# W-method solvers and strict Rosenbrock solvers
# ============================================================================

w_methods = [
    Rosenbrock23(),
    Rosenbrock32(),
    Rodas23W(),
    ROS2S(),
    ROS34PW1a(),
    ROS34PW2(),
    ROS34PW3(),
    ROK4a(),
]

strict_rosenbrock = [
    ROS3P(),
    Rodas3(),
    Rodas3P(),
    Rodas4(),
    Rodas4P(),
    Rodas4P2(),
    Rodas5(),
    Rodas5P(),
    Rodas5Pe(),
    Rodas5Pr(),
]

@testset "Jacobian Reuse for Rosenbrock-W Methods" begin

    # ========================================================================
    @testset "isWmethod trait consistency" begin
        # Verify W-methods have the trait set
        for alg in w_methods
            @test OrdinaryDiffEqRosenbrock.isWmethod(alg) == true
        end
        # Verify strict Rosenbrock methods don't have it
        for alg in strict_rosenbrock
            @test OrdinaryDiffEqRosenbrock.isWmethod(alg) == false
        end
    end

    # ========================================================================
    @testset "JacReuseState construction" begin
        jr = JacReuseState(0.0)
        @test jr.last_dtgamma == 0.0
        @test jr.pending_dtgamma == 0.0
        @test jr.last_naccept == 0
        @test jr.max_jac_age == 50
        @test jr.cached_J === nothing
        @test jr.cached_dT === nothing
    end

    # ========================================================================
    @testset "Convergence preserved with J reuse - W-methods" begin
        dts = (1 / 2) .^ (6:-1:3)
        testTol = 0.5

        # Rosenbrock23 is order 2
        order = estimate_order(linear2d_prob, Rosenbrock23(), dts, exact_2d)
        @test order ≈ 2 atol = testTol

        # Rosenbrock32 converges at order 3
        order = estimate_order(linear2d_prob, Rosenbrock32(), dts, exact_2d)
        @test order ≈ 3 atol = testTol

        # Rodas23W is 2/3 pair; fixed-step convergence gives order 2
        order = estimate_order(linear2d_prob, Rodas23W(), dts, exact_2d)
        @test order ≈ 2 atol = testTol

        # ROS34PW3 is order 4
        order = estimate_order(linear2d_prob, ROS34PW3(), dts, exact_2d)
        @test order ≈ 4 atol = testTol
    end

    # ========================================================================
    @testset "Convergence preserved - strict Rosenbrock" begin
        dts = (1 / 2) .^ (6:-1:3)
        testTol = 0.5

        for alg in [ROS3P(), Rodas3(), Rodas4(), Rodas5()]
            order_expected = OrdinaryDiffEqRosenbrock.alg_order(alg)
            order = estimate_order(linear2d_prob, alg, dts, exact_2d)
            @test order ≈ order_expected atol = testTol
        end
    end

    # ========================================================================
    @testset "Jacobian count reduced for W-methods on stiff problems" begin
        # W-methods should reuse Jacobians, resulting in njacs < naccept
        for alg in [Rosenbrock23(), Rodas23W(), ROS34PW3()]
            sol = solve(vdp_prob, alg, reltol = 1.0e-6, abstol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            # W-methods should have fewer Jacobian evaluations than accepted steps
            if sol.stats.naccept > 10
                @test sol.stats.njacs < sol.stats.naccept
            end
        end
    end

    # ========================================================================
    @testset "Strict Rosenbrock methods compute J every step" begin
        # Strict Rosenbrock methods should compute J on every (non-repeated) step
        for alg in [Rodas3(), Rodas4(), Rodas5(), Rodas5P()]
            sol = solve(vdp_prob, alg, reltol = 1.0e-6, abstol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            # For strict Rosenbrock, njacs should be >= naccept
            if sol.stats.naccept > 5
                @test sol.stats.njacs >= sol.stats.naccept
            end
        end
    end

    # ========================================================================
    @testset "ROBER benchmark accuracy" begin
        # Solve ROBER with reference method
        ref_sol = solve(rober_prob, Rodas5P(), reltol = 1.0e-12, abstol = 1.0e-12)

        for alg in [Rosenbrock23(), Rodas23W(), ROS34PW3()]
            sol = solve(rober_prob, alg, reltol = 1.0e-6, abstol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            # Check solution is reasonably close at the endpoint
            @test norm(sol.u[end] - ref_sol.u[end]) / norm(ref_sol.u[end]) < 1.0e-3
        end
    end

    # ========================================================================
    @testset "Van der Pol benchmark accuracy and savings" begin
        ref_sol = solve(vdp_prob, Rodas5P(), reltol = 1.0e-12, abstol = 1.0e-12)

        for alg in [Rosenbrock23(), Rodas23W()]
            sol = solve(vdp_prob, alg, reltol = 1.0e-6, abstol = 1.0e-8)
            @test SciMLBase.successful_retcode(sol)
            @test norm(sol.u[end] - ref_sol.u[end]) / norm(ref_sol.u[end]) < 1.0e-2
        end
    end

    # ========================================================================
    @testset "Mass matrix DAE with W-method" begin
        # x' = -x + y, 0 = x - y  =>  y = x, x' = 0  =>  x(t) = y(t) = 1
        for alg in [Rodas23W(), ROS34PW3()]
            sol = solve(
                mass_prob, alg, reltol = 1.0e-8, abstol = 1.0e-10,
                initializealg = SciMLBase.CheckInit()
            )
            @test SciMLBase.successful_retcode(sol)
            # The algebraic constraint x = y should be satisfied
            @test abs(sol.u[end][1] - sol.u[end][2]) < 1.0e-5
            # x(t) = 1 (constant solution)
            @test abs(sol.u[end][1] - 1.0) < 1.0e-4
        end
    end

    # ========================================================================
    @testset "Exponential decay - all solvers produce correct result" begin
        for alg in [w_methods; strict_rosenbrock]
            sol = solve(expdecay_prob, alg, reltol = 1.0e-8, abstol = 1.0e-10)
            @test SciMLBase.successful_retcode(sol)
            # Exact solution: u(t) = u0 * exp(-t)
            exact = [exp(-1.0), 2 * exp(-1.0), 3 * exp(-1.0)]
            # Use looser tolerance to accommodate low-order methods (order 2)
            @test norm(sol.u[end] - exact) < 1.0e-2
        end
    end
end

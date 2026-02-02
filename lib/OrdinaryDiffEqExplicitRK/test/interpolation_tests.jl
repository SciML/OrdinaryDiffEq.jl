using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, constructDormandPrince
using OrdinaryDiffEqCore
using DiffEqBase
using Test
import SciMLBase

# ============================================================================
# Test Problems
# ============================================================================

f_linear(u, p, t) = 1.01 * u
prob_ode_linear = ODEProblem(
    ODEFunction(f_linear; analytic = (u0, p, t) -> u0 * exp(1.01t)),
    1 / 2, (0.0, 1.0)
)

function f_2Dlinear!(du, u, p, t)
    du[1] = 1.01 * u[1]
    return du[2] = 1.01 * u[2]
end
prob_ode_2Dlinear = ODEProblem(
    ODEFunction(f_2Dlinear!; analytic = (u0, p, t) -> u0 .* exp(1.01t)),
    [1 / 2, 1 / 2], (0.0, 1.0)
)

# ============================================================================
# Basic Solve Tests
# ============================================================================

@testset "Basic Solve" begin
    for tableau_fn in [constructTsit5ExplicitRK, constructDormandPrince]
        sol = solve(prob_ode_linear, ExplicitRK(tableau = tableau_fn()))
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)

        sol = solve(prob_ode_2Dlinear, ExplicitRK(tableau = tableau_fn()))
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)
    end
end

# ============================================================================
# Interpolation Tests
# ============================================================================

@testset "Interpolation" begin
    f_decay(u, p, t) = -u
    prob_interp = ODEProblem(
        ODEFunction(f_decay; analytic = (u0, p, t) -> exp(-t)), 1.0, (0.0, 1.0)
    )

    sol = solve(prob_interp, ExplicitRK(tableau = constructTsit5ExplicitRK()), dense = true)
    @test SciMLBase.successful_retcode(sol)
    @test sol(0.0) ≈ 1.0 atol = 1.0e-8
    @test sol(0.5) ≈ exp(-0.5) atol = 0.1
    @test sol(1.0) ≈ exp(-1.0) atol = 0.1

    # Vector problem (in-place)
    function f_vec!(du, u, p, t)
        du[1] = -u[1]
        du[2] = -2 * u[2]
    end
    exact_vec(t) = [exp(-t), exp(-2t)]
    prob_vec = ODEProblem(
        ODEFunction(f_vec!; analytic = (u0, p, t) -> exact_vec(t)),
        [1.0, 1.0], (0.0, 1.0)
    )

    sol = solve(prob_vec, ExplicitRK(tableau = constructTsit5ExplicitRK()), dense = true)
    @test SciMLBase.successful_retcode(sol)
    @test sol(0.5) ≈ exact_vec(0.5) atol = 0.1

    # Interpolation with idxs
    @test sol(0.5, idxs = 1) ≈ exp(-0.5) atol = 0.1
    @test sol(0.5, idxs = 2) ≈ exp(-1.0) atol = 0.1
    @test sol(0.5, idxs = [1, 2]) ≈ exact_vec(0.5) atol = 0.1
end

# ============================================================================
# L2 Convergence Order Tests
# ============================================================================

function compute_midstep_error(sol, exact_fn)
    max_err = 0.0
    for i in 1:(length(sol.t) - 1)
        t_mid = (sol.t[i] + sol.t[i + 1]) / 2
        interp_val = sol(t_mid)
        exact_val = exact_fn(t_mid)
        if interp_val isa Number
            err = abs(interp_val - exact_val)
        else
            err = maximum(abs.(interp_val .- exact_val))
        end
        max_err = max(max_err, err)
    end
    return max_err
end

function estimate_order(errors, dts)
    orders = Float64[]
    for i in 2:length(errors)
        order = log(errors[i - 1] / errors[i]) / log(dts[i - 1] / dts[i])
        push!(orders, order)
    end
    return orders
end

@testset "L2 Convergence - Scalar" begin
    f_conv(u, p, t) = -u
    exact_scalar(t) = exp(-t)
    prob_conv = ODEProblem(
        ODEFunction(f_conv; analytic = (u0, p, t) -> exact_scalar(t)), 1.0, (0.0, 1.0)
    )

    tableau = constructTsit5ExplicitRK()
    dts = [1 / 2^k for k in 2:6]
    errors = Float64[]

    for dt in dts
        sol = solve(
            prob_conv, ExplicitRK(tableau = tableau); dt = dt, adaptive = false, dense = true
        )
        push!(errors, compute_midstep_error(sol, exact_scalar))
    end

    orders = estimate_order(errors, dts)
    avg_order = sum(orders) / length(orders)
    @test avg_order > 3.5
end

@testset "L2 Convergence - Vector" begin
    function f_conv_vec!(du, u, p, t)
        du[1] = -u[1]
        du[2] = -2 * u[2]
    end
    exact_vector(t) = [exp(-t), exp(-2t)]
    prob_conv_vec = ODEProblem(
        ODEFunction(f_conv_vec!; analytic = (u0, p, t) -> exact_vector(t)),
        [1.0, 1.0], (0.0, 1.0)
    )

    tableau = constructTsit5ExplicitRK()
    dts = [1 / 2^k for k in 2:6]
    errors = Float64[]

    for dt in dts
        sol = solve(
            prob_conv_vec, ExplicitRK(tableau = tableau); dt = dt, adaptive = false, dense = true
        )
        err = compute_midstep_error(sol, exact_vector)
        push!(errors, err)
    end

    orders = estimate_order(errors, dts)
    avg_order = sum(orders) / length(orders)
    @test avg_order > 3.5
end

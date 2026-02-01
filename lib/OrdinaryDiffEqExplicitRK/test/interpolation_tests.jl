"""
    ExplicitRK Tests -

"""

using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, constructDormandPrince
using OrdinaryDiffEqCore
using DiffEqBase
using Test
using Random
import SciMLBase

Random.seed!(100)

# ============================================================================
# Test Problems
# ============================================================================

# Scalar linear ODE: y' = 1.01y
f_linear(u, p, t) = 1.01 * u
prob_ode_linear = ODEProblem(
    ODEFunction(f_linear; analytic=(u0, p, t) -> u0 * exp(1.01t)),
    1/2, (0.0, 1.0)
)

# 2D linear ODE (in-place)
function f_2Dlinear!(du, u, p, t)
    du[1] = 1.01 * u[1]
    du[2] = 1.01 * u[2]
end
prob_ode_2Dlinear = ODEProblem(
    ODEFunction(f_2Dlinear!; analytic=(u0, p, t) -> u0 .* exp(1.01t)),
    [1/2, 1/2], (0.0, 1.0)
)

# ============================================================================
# Basic Solve Tests
# ============================================================================

prob = prob_ode_linear
sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRK()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

sol = solve(prob, ExplicitRK(tableau=constructDormandPrince()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

prob = prob_ode_2Dlinear
sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRK()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

sol = solve(prob, ExplicitRK(tableau=constructDormandPrince()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

# ============================================================================
# Tableau Comparison Tests
# ============================================================================

# Compare Tsit5 vs Dormand-Prince tableau
prob = prob_ode_linear
tableau_tsit5 = constructTsit5ExplicitRK()
tableau_dp = constructDormandPrince()

sol1 = solve(prob, ExplicitRK(tableau=tableau_tsit5); dt=1/2^6, adaptive=false, save_everystep=false)
sol2 = solve(prob, ExplicitRK(tableau=tableau_dp); dt=1/2^6, adaptive=false, save_everystep=false)
# Both are 5th order methods, should give similar (but not identical) results
@test abs(sol1.u[end] - sol2.u[end]) < 1e-4

prob = prob_ode_2Dlinear
sol1 = solve(prob, ExplicitRK(tableau=tableau_tsit5); dt=1/2^3, adaptive=false, save_everystep=false)
sol2 = solve(prob, ExplicitRK(tableau=tableau_dp); dt=1/2^3, adaptive=false, save_everystep=false)
@test maximum(abs.(sol1.u[end] - sol2.u[end])) < 1e-3

# Adaptive stepping should work for both
sol1 = solve(prob, ExplicitRK(tableau=tableau_tsit5); dt=1/2^6)
sol2 = solve(prob, ExplicitRK(tableau=tableau_dp); dt=1/2^6)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# ============================================================================
# Interpolation Tests
# ============================================================================

# Test problem with known analytical solution
f(u, p, t) = -u
prob_interp = ODEProblem(ODEFunction(f; analytic=(u0, p, t) -> exp(-t)), 1.0, (0.0, 1.0))


sol = solve(prob_interp, ExplicitRK(tableau=constructTsit5ExplicitRK()), dense=true)
@test SciMLBase.successful_retcode(sol)
@test sol(0.0) ≈ 1.0 atol=1e-8
@test sol(0.5) ≈ exp(-0.5) atol=0.1
@test sol(1.0) ≈ exp(-1.0) atol=0.1


function f_vec!(du, u, p, t)
    du[1] = -u[1]
    du[2] = -2*u[2]
end
exact_vec(t) = [exp(-t), exp(-2t)]
prob_interp_vec = ODEProblem(ODEFunction(f_vec!; analytic=(u0, p, t) -> exact_vec(t)), [1.0, 1.0], (0.0, 1.0))

sol = solve(prob_interp_vec, ExplicitRK(tableau=constructTsit5ExplicitRK()), dense=true)
@test SciMLBase.successful_retcode(sol)
@test sol(0.5) ≈ exact_vec(0.5) atol=0.1

# Interpolation with idxs
@test sol(0.5, idxs=1) ≈ exp(-0.5) atol=0.1
@test sol(0.5, idxs=2) ≈ exp(-1.0) atol=0.1
@test sol(0.5, idxs=[1,2]) ≈ exact_vec(0.5) atol=0.1

# Interpolation at step boundaries
@test sol(0.0) ≈ [1.0, 1.0] atol=1e-6
@test sol(1.0) ≈ sol.u[end] atol=1e-5

# ============================================================================
# L2 Convergence Order Tests
# ============================================================================

"""
Compute the L2 error of the interpolant over the interval [t0, t1].

"""
function compute_L2_error(sol, exact_fn; npoints=1000)
    t0, t1 = sol.prob.tspan
    ts = range(t0, t1, length=npoints)
    dt_spacing = (t1 - t0) / (npoints - 1)

    # Compute L2 error using trapezoidal-like approximation
    total_sq_error = 0.0
    for t in ts
        interp_val = sol(t)
        exact_val = exact_fn(t)
        if interp_val isa Number
            total_sq_error += (interp_val - exact_val)^2
        else
            total_sq_error += sum((interp_val .- exact_val).^2)
        end
    end
    return sqrt(total_sq_error * dt_spacing)
end

"""
Compute max error at intermediate points (mid-step) only.
"""
function compute_midstep_error(sol, exact_fn)
    max_err = 0.0
    for i in 1:(length(sol.t)-1)
        t_mid = (sol.t[i] + sol.t[i+1]) / 2
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

"""
Estimate the convergence order from errors at successive step sizes.
"""
function estimate_order(errors, dts)
    orders = Float64[]
    for i in 2:length(errors)
        order = log(errors[i-1] / errors[i]) / log(dts[i-1] / dts[i])
        push!(orders, order)
    end
    return orders
end

# Test L2 convergence for scalar problem
@testset "L2 Convergence - Scalar" begin
    f_conv(u, p, t) = -u
    exact_scalar(t) = exp(-t)
    prob_conv = ODEProblem(ODEFunction(f_conv; analytic=(u0, p, t) -> exact_scalar(t)), 1.0, (0.0, 1.0))

    tableau = constructTsit5ExplicitRK()
    dts = [1/2^k for k in 2:6]  # h = 1/4, 1/8, 1/16, 1/32, 1/64
    errors_midstep = Float64[]

    for dt in dts
        sol = solve(prob_conv, ExplicitRK(tableau=tableau); dt=dt, adaptive=false, dense=true)
        push!(errors_midstep, compute_midstep_error(sol, exact_scalar))
    end

    orders_midstep = estimate_order(errors_midstep, dts)

    # Tsit5 interpolation should achieve at least 4th order convergence
    avg_order = sum(orders_midstep) / length(orders_midstep)
    @test avg_order > 3.5
end

# Test L2 convergence for vector problem (in-place)
@testset "L2 Convergence - Vector" begin
    function f_conv_vec!(du, u, p, t)
        du[1] = -u[1]
        du[2] = -2*u[2]
    end
    exact_vector(t) = [exp(-t), exp(-2t)]
    prob_conv_vec = ODEProblem(ODEFunction(f_conv_vec!; analytic=(u0, p, t) -> exact_vector(t)),
                               [1.0, 1.0], (0.0, 1.0))

    tableau = constructTsit5ExplicitRK()
    dts = [1/2^k for k in 2:6]
    errors = Float64[]

    for dt in dts
        sol = solve(prob_conv_vec, ExplicitRK(tableau=tableau); dt=dt, adaptive=false, dense=true)
        err = compute_L2_error(sol, exact_vector)
        push!(errors, err)
    end

    orders = estimate_order(errors, dts)
    avg_order = sum(orders) / length(orders)
    @test avg_order > 3.5
end

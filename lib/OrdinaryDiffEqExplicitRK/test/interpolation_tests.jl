using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, constructDormandPrince,
    ExplicitRKConstantCache
using OrdinaryDiffEqCore
using DiffEqBase
using ForwardDiff
using Test
import SciMLBase
import JLArrays

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
    du[2] = 1.01 * u[2]
    return
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
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)

        sol = solve(prob_ode_2Dlinear, ExplicitRK(tableau = tableau_fn()))
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end
end

@testset "Float32 in-place ExplicitRK" begin
    function f32!(du, u, p, t)
        du[1] = 1.01f0 * u[1]
        return nothing
    end
    prob = ODEProblem(f32!, Float32[0.5], (0.0f0, 1.0f0))
    sol = solve(prob, ExplicitRK())
    @test SciMLBase.successful_retcode(sol)
    @test eltype(sol.u[end]) === Float32
    @test eltype(sol.t) === Float32
end

@testset "ExplicitRKConstantCache c eltype inference" begin
    tab = constructDormandPrince(BigFloat)
    cache = @inferred ExplicitRKConstantCache(tab, [1.0], Float64)
    @test eltype(cache.c) === BigFloat
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
            prob_conv, ExplicitRK(; tableau); dt, adaptive = false, dense = true
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
            prob_conv_vec, ExplicitRK(; tableau); dt, adaptive = false, dense = true
        )
        err = compute_midstep_error(sol, exact_vector)
        push!(errors, err)
    end

    orders = estimate_order(errors, dts)
    avg_order = sum(orders) / length(orders)
    @test avg_order > 3.5
end

@testset "interpolating before the first step" begin
    f!(du, u, p, t) = (du .= -u; nothing)
    prob = ODEProblem(f!, [1.0, 2.0], (0.0, 1.0))

    integ = init(prob, ExplicitRK(); dt = 0.1)
    @test integ.sol(0.0) ≈ [1.0, 2.0]

    integ.sol(0.0, Val{1})
    @test integ.sol.k[1][1] ≈ [-1.0, -2.0]
    @test integ.sol.k[1][2] ≈ [-1.0, -2.0]

    integ2 = init(prob, ExplicitRK(); dt = 0.1)
    out = zeros(2)
    integ2.sol(out, 0.0, Val{1})
    @test integ2.sol.k[1][1] ≈ [-1.0, -2.0]
    @test integ2.sol.k[1][2] ≈ [-1.0, -2.0]
end

@testset "fallback compute_stages! (>17 stages) on JLArray" begin
    n = 19
    A = zeros(n, n)
    for i in 2:n, j in 1:(i - 1)
        A[i, j] = 1 / (i - 1)
    end
    c = [0; fill(0.5, n - 1)]
    α = zeros(n)
    α[end] = 1.0
    αEEst = zeros(n)
    αEEst[end] = 0.5
    αEEst[1] = -0.5
    tab = DiffEqBase.ExplicitRKTableau(A, c, α, 2; αEEst, adaptiveorder = 1)
    alg = ExplicitRK(tableau = tab)
    f_jl!(du, u, p, t) = (du .= -0.5 .* u; nothing)

    sol_cpu = solve(
        ODEProblem(f_jl!, ones(10), (0.0, 1.0)), alg;
        adaptive = false, dt = 0.1, dense = false,
    )
    @test SciMLBase.successful_retcode(sol_cpu)

    JLArrays.allowscalar(false)
    try
        sol_jl = solve(
            ODEProblem(f_jl!, JLArrays.JLVector(ones(10)), (0.0, 1.0)), alg;
            adaptive = false, dt = 0.1, dense = false,
        )
        @test SciMLBase.successful_retcode(sol_jl)
        @test Array(sol_jl.u[end]) == sol_cpu.u[end]
    finally
        JLArrays.allowscalar(true)
    end
end

@testset "idxs interpolation under ForwardDiff Dual times" begin
    function f_vec!(du, u, p, t)
        du[1] = -u[1]
        du[2] = -2 * u[2]
        return nothing
    end
    prob = ODEProblem(f_vec!, [1.0, 1.0], (0.0, 1.0))
    sol = solve(prob, ExplicitRK(tableau = constructTsit5ExplicitRK()); dense = true)
    @test SciMLBase.successful_retcode(sol)

    t0 = sol.t[3]
    expected_d1 = [-exp(-t0), -2 * exp(-2t0)]
    # Reference for the derivative-of-derivative check: idxs must select the
    # same components as the full-state interpolant derivative.
    expected_d2 = ForwardDiff.derivative(t -> sol(t, Val{1}), t0)

    # sol(t; idxs) with an index vector routes through the in-place kernel
    # internally, so it exercises the same code path as sol(out, t; idxs).
    @test ForwardDiff.derivative(t -> sol(t; idxs = 1), t0) ≈ expected_d1[1] atol = 1.0e-6
    @test ForwardDiff.derivative(t -> sol(t; idxs = [1, 2]), t0) ≈ expected_d1 atol = 1.0e-6
    @test ForwardDiff.derivative(t -> sol(t, Val{1}; idxs = [1, 2]), t0) ≈ expected_d2

    @test ForwardDiff.derivative(t0) do t
        out = zeros(typeof(t), 2)
        sol(out, t; idxs = [1, 2])
        out
    end ≈ expected_d1 atol = 1.0e-6

    @test ForwardDiff.derivative(t0) do t
        out = zeros(typeof(t), 1)
        sol(out, t; idxs = 1)
        out[1]
    end ≈ expected_d1[1] atol = 1.0e-6

    @test ForwardDiff.derivative(t0) do t
        out = zeros(typeof(t), 2)
        sol(out, t, Val{1}; idxs = [1, 2])
        out
    end ≈ expected_d2
end

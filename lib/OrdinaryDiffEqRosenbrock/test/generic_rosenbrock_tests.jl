using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqRosenbrock: constructRodas4, constructRodas42, constructRodas4P,
    constructRodas4P2, constructRodas5, constructRodas5P, constructRodas6P
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
# Helper Functions
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

# ============================================================================
# Tableau Metadata for Testing
# ============================================================================

const TABLEAU_CONFIGS = [
    (name = "Rodas4", constructor = constructRodas4, order = 4, sol_stages = 6),
    (name = "Rodas42", constructor = constructRodas42, order = 4, sol_stages = 6),
    (name = "Rodas4P", constructor = constructRodas4P, order = 4, sol_stages = 6),
    (name = "Rodas4P2", constructor = constructRodas4P2, order = 4, sol_stages = 6),
    (name = "Rodas5", constructor = constructRodas5, order = 5, sol_stages = 8),
    (name = "Rodas5P", constructor = constructRodas5P, order = 5, sol_stages = 8),
    (name = "Rodas6P", constructor = constructRodas6P, order = 6, sol_stages = 16),
]

function make_generic(config; kwargs...)
    GenericRosenbrock(;
        tableau = config.constructor(), order = config.order,
        num_solution_stages = config.sol_stages, kwargs...
    )
end
function make_generic(config, ::Type{T}; kwargs...) where {T}
    GenericRosenbrock(;
        tableau = config.constructor(T), order = config.order,
        num_solution_stages = config.sol_stages, kwargs...
    )
end

# ============================================================================
# GenericRosenbrock Basic Solve Tests
# ============================================================================

@testset "GenericRosenbrock Basic Solve" begin
    @testset "$(config.name)" for config in TABLEAU_CONFIGS
        @testset "Scalar problem" begin
            sol = solve(prob_ode_linear, make_generic(config))
            @test SciMLBase.successful_retcode(sol)
        end

        @testset "Vector problem" begin
            sol = solve(prob_ode_2Dlinear, make_generic(config))
            @test SciMLBase.successful_retcode(sol)
        end
    end
end

# ============================================================================
# GenericRosenbrock Interpolation Tests
# ============================================================================

@testset "GenericRosenbrock Interpolation" begin
    @testset "$(config.name)" for config in TABLEAU_CONFIGS
        f_decay(u, p, t) = -u
        prob_interp = ODEProblem(
            ODEFunction(f_decay; analytic = (u0, p, t) -> exp(-t)), 1.0, (0.0, 1.0)
        )

        @testset "Scalar interpolation" begin
            sol = solve(
                prob_interp, make_generic(config), dense = true
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol(0.0) ≈ 1.0 atol = 1.0e-8
            @test sol(0.5) ≈ exp(-0.5) atol = 0.01
            @test sol(1.0) ≈ exp(-1.0) atol = 0.01
        end

        @testset "Vector interpolation" begin
            function f_vec!(du, u, p, t)
                du[1] = -u[1]
                du[2] = -2 * u[2]
            end
            exact_vec(t) = [exp(-t), exp(-2t)]
            prob_vec = ODEProblem(
                ODEFunction(f_vec!; analytic = (u0, p, t) -> exact_vec(t)),
                [1.0, 1.0], (0.0, 1.0)
            )

            sol = solve(
                prob_vec, make_generic(config), dense = true
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol(0.5) ≈ exact_vec(0.5) atol = 0.01
        end

        @testset "Interpolation with idxs" begin
            function f_vec!(du, u, p, t)
                du[1] = -u[1]
                du[2] = -2 * u[2]
            end
            exact_vec(t) = [exp(-t), exp(-2t)]
            prob_vec = ODEProblem(
                ODEFunction(f_vec!; analytic = (u0, p, t) -> exact_vec(t)),
                [1.0, 1.0], (0.0, 1.0)
            )

            sol = solve(
                prob_vec, make_generic(config), dense = true
            )
            @test sol(0.5, idxs = 1) ≈ exp(-0.5) atol = 0.01
            @test sol(0.5, idxs = 2) ≈ exp(-1.0) atol = 0.01
            @test sol(0.5, idxs = [1, 2]) ≈ exact_vec(0.5) atol = 0.01
        end
    end
end

# ============================================================================
# GenericRosenbrock L2 Convergence Tests
# ============================================================================

@testset "GenericRosenbrock L2 Convergence - Scalar" begin
    @testset "$(config.name)" for config in TABLEAU_CONFIGS
        f_conv(u, p, t) = -u
        exact_scalar(t) = exp(-t)
        prob_conv = ODEProblem(
            ODEFunction(f_conv; analytic = (u0, p, t) -> exact_scalar(t)), 1.0, (0.0, 1.0)
        )

        dts = [1 / 2^k for k in 2:6]
        errors = Float64[]

        for dt in dts
            sol = solve(
                prob_conv, make_generic(config);
                dt = dt, adaptive = false, dense = true
            )
            push!(errors, compute_midstep_error(sol, exact_scalar))
        end

        orders = estimate_order(errors, dts)
        avg_order = sum(orders) / length(orders)
        expected_order = config.order - 0.5
        @test avg_order > expected_order
    end
end

@testset "GenericRosenbrock L2 Convergence - Vector" begin
    @testset "$(config.name)" for config in TABLEAU_CONFIGS
        function f_conv_vec!(du, u, p, t)
            du[1] = -u[1]
            du[2] = -2 * u[2]
        end
        exact_vector(t) = [exp(-t), exp(-2t)]
        prob_conv_vec = ODEProblem(
            ODEFunction(f_conv_vec!; analytic = (u0, p, t) -> exact_vector(t)),
            [1.0, 1.0], (0.0, 1.0)
        )

        dts = [1 / 2^k for k in 2:6]
        errors = Float64[]

        for dt in dts
            sol = solve(
                prob_conv_vec, make_generic(config);
                dt = dt, adaptive = false, dense = true
            )
            err = compute_midstep_error(sol, exact_vector)
            push!(errors, err)
        end

        orders = estimate_order(errors, dts)
        avg_order = sum(orders) / length(orders)
        expected_order = config.order - 0.5
        @test avg_order > expected_order
    end
end

# ============================================================================
# GenericRosenbrock vs Specialized Solver Comparison
# ============================================================================

@testset "GenericRosenbrock Matches Specialized Solvers" begin
    f_test(u, p, t) = -u
    prob = ODEProblem(f_test, 1.0, (0.0, 2.0))
    test_times = [0.0, 0.5, 1.0, 1.5, 2.0]

    @testset "Rodas4" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas4()))
        sol_specialized = solve(prob, Rodas4())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end

    @testset "Rodas42" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas42()))
        sol_specialized = solve(prob, Rodas42())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-5
        end
    end

    @testset "Rodas4P" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas4P()))
        sol_specialized = solve(prob, Rodas4P())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end

    @testset "Rodas4P2" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas4P2()))
        sol_specialized = solve(prob, Rodas4P2())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end

    @testset "Rodas5" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas5()))
        sol_specialized = solve(prob, Rodas5())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end

    @testset "Rodas5P" begin
        sol_generic = solve(prob, GenericRosenbrock(tableau = constructRodas5P()))
        sol_specialized = solve(prob, Rodas5P())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end

    @testset "Rodas6P" begin
        sol_generic = solve(prob, GenericRosenbrock(
            tableau = constructRodas6P(), order = 6, num_solution_stages = 16
        ))
        sol_specialized = solve(prob, Rodas6P())

        @test SciMLBase.successful_retcode(sol_generic)
        @test SciMLBase.successful_retcode(sol_specialized)

        for t in test_times
            @test sol_generic(t) ≈ sol_specialized(t) atol = 1e-6
        end
    end
end

# ============================================================================
# GenericRosenbrock Dense Output Tests
# ============================================================================

@testset "GenericRosenbrock Dense Output" begin
    @testset "$(config.name)" for config in TABLEAU_CONFIGS
        f_decay(u, p, t) = -u
        prob = ODEProblem(f_decay, 1.0, (0.0, 1.0))

        sol = solve(prob, make_generic(config), dense = true)

        @testset "Interpolation at multiple points" begin
            test_times = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
            for t in test_times
                @test sol(t) ≈ exp(-t) atol = 1e-4
            end
        end

        @testset "Midstep interpolation" begin
            for i in 1:(length(sol.t) - 1)
                t_mid = (sol.t[i] + sol.t[i + 1]) / 2
                @test sol(t_mid) ≈ exp(-t_mid) atol = 1e-3
            end
        end
    end
end

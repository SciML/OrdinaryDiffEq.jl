using DelayDiffEq
using OrdinaryDiffEqTsit5
using Test

# Test allocation regression for key operations
# These tests verify that in-place variants of critical operations
# remain allocation-free (or have minimal allocations)

@testset "Allocation Regression Tests" begin
    # Setup a simple DDE problem
    function f!(du, u, h, p, t)
        du[1] = -h(p, t - 1.0)[1]
    end

    function h_inplace!(val, p, t; idxs = nothing)
        if idxs === nothing
            val[1] = 1.0
        else
            val[idxs] = 1.0
        end
        return val
    end

    function h(p, t; idxs = nothing)
        if idxs === nothing
            return [1.0]
        else
            return 1.0
        end
    end

    u0 = [1.0]
    tspan = (0.0, 5.0)
    p = nothing

    prob = DDEProblem(f!, u0, h, tspan, p; constant_lags = [1.0])
    alg = MethodOfSteps(Tsit5())

    # Solve to have a dense history
    sol = solve(prob, alg)

    @testset "Solution interpolation (in-place)" begin
        # In-place interpolation should be allocation-free after warmup
        u_cache = zeros(1)

        # Warm up
        sol(u_cache, 2.5)

        # Test for zero allocations (or very minimal)
        allocs = @allocated sol(u_cache, 2.5)
        @test allocs == 0
    end

    @testset "History function (in-place)" begin
        # Create integrator to access history function
        integrator = init(prob, alg)

        hf = integrator.history
        cache_val = zeros(1)

        # Warm up
        hf(cache_val, p, 0.5)

        # In-place history function should be allocation-free
        allocs = @allocated hf(cache_val, p, 0.5)
        @test allocs == 0
    end

    @testset "Solve allocation bounds" begin
        # Test that total solve allocations stay within reasonable bounds
        # This helps detect allocation regressions in the solve path

        # Warm up
        solve(prob, alg)

        # Count allocations for a simple problem
        allocs = @allocated solve(prob, alg)

        # The solve should allocate less than 100KB for this simple problem
        # This is a regression test - if allocations increase significantly,
        # the test will fail
        @test allocs < 100_000

        # Also check number of allocations stays bounded
        # (this is more stable than bytes which depend on array sizes)
        num_allocs = 0
        solve(prob, alg)  # warmup
        stats = @timed solve(prob, alg)
        # We don't have direct access to allocation count in @timed,
        # but we can verify the bytes are bounded
        @test stats.bytes < 100_000
    end

    @testset "Step execution allocation stability" begin
        # For a simple DDE, each step should have bounded allocations
        integrator = init(prob, alg)

        # Take a few steps to warm up
        for _ in 1:5
            step!(integrator)
        end

        # Reset and measure
        reinit!(integrator)
        step!(integrator)  # First step after reinit

        # Measure subsequent step allocations
        allocs_per_step = @allocated step!(integrator)

        # Each step should allocate less than 5KB
        # (most allocations are for growing solution arrays)
        @test allocs_per_step < 5_000
    end

    @testset "Multi-dimensional problem allocations" begin
        # Test allocation bounds for multi-dimensional problems
        function f_multi!(du, u, h, p, t)
            h_val = h(p, t - 1.0)
            for i in 1:6
                du[i] = -0.1 * u[i] + 0.5 * h_val[i]
            end
        end

        function h_multi(p, t; idxs = nothing)
            if idxs === nothing
                return ones(6)
            else
                return 1.0
            end
        end

        u0_multi = ones(6)
        prob_multi = DDEProblem(f_multi!, u0_multi, h_multi, tspan, p; constant_lags = [1.0])

        # Warm up
        sol_multi = solve(prob_multi, alg)
        solve(prob_multi, alg)

        # Multi-dimensional solve should also have bounded allocations
        allocs = @allocated solve(prob_multi, alg)
        # Allow more allocations for larger state vectors, but still bounded
        @test allocs < 150_000

        # Test in-place interpolation for multi-dimensional
        u_cache_multi = zeros(6)
        sol_multi(u_cache_multi, 2.5)
        allocs_interp = @allocated sol_multi(u_cache_multi, 2.5)
        @test allocs_interp == 0
    end

    @testset "Scalar (out-of-place) problem allocations" begin
        # Test allocation bounds for scalar out-of-place problems
        function f_scalar(u, h, p, t)
            return -h(p, t - 1.0)
        end

        function h_scalar(p, t; idxs = nothing)
            return 1.0
        end

        prob_scalar = DDEProblem(f_scalar, 1.0, h_scalar, tspan, p; constant_lags = [1.0])

        # Warm up
        sol_scalar = solve(prob_scalar, alg)
        solve(prob_scalar, alg)

        # Scalar problem should have lower allocations
        allocs = @allocated solve(prob_scalar, alg)
        @test allocs < 50_000

        # Test step allocations for scalar
        integrator_scalar = init(prob_scalar, alg)
        for _ in 1:5
            step!(integrator_scalar)
        end
        reinit!(integrator_scalar)
        step!(integrator_scalar)

        allocs_step = @allocated step!(integrator_scalar)
        # Scalar steps should allocate less than vector steps
        @test allocs_step < 1_000
    end

    @testset "Dependent delay allocations" begin
        # Test allocation bounds for problems with dependent delays
        function f_dep!(du, u, h, p, t)
            tau = 0.5 + 0.1 * abs(u[1])
            h_val = h(p, t - tau)
            du[1] = -h_val[1]
        end

        function h_dep(p, t; idxs = nothing)
            if idxs === nothing
                return [1.0]
            else
                return 1.0
            end
        end

        prob_dep = DDEProblem(
            f_dep!, [1.0], h_dep, (0.0, 3.0), nothing;
            dependent_lags = [(u, p, t) -> 0.5 + 0.1 * abs(u[1])]
        )

        # Warm up
        sol_dep = solve(prob_dep, alg)
        solve(prob_dep, alg)

        # Dependent delays have more overhead but should still be bounded
        allocs = @allocated solve(prob_dep, alg)
        @test allocs < 100_000
    end

    @testset "Unconstrained mode allocations" begin
        # Test allocation bounds for unconstrained mode (fixed-point iteration)
        alg_uc = MethodOfSteps(Tsit5(); constrained = false)

        # Warm up
        sol_uc = solve(prob, alg_uc)
        solve(prob, alg_uc)

        # Unconstrained mode should have similar allocation bounds
        allocs = @allocated solve(prob, alg_uc)
        @test allocs < 100_000

        # Test step allocations in unconstrained mode
        integrator_uc = init(prob, alg_uc)
        for _ in 1:5
            step!(integrator_uc)
        end
        reinit!(integrator_uc)
        step!(integrator_uc)

        allocs_step = @allocated step!(integrator_uc)
        @test allocs_step < 5_000
    end
end

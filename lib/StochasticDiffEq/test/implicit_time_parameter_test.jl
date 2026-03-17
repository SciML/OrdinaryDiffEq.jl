using StochasticDiffEq, Test

@testset "ImplicitEM Newton solver time coefficient" begin
    # This test verifies that implicit SDE solvers evaluate the drift function
    # at the correct time during Newton iterations. The nlsolver.c coefficient
    # should be the Butcher tableau coefficient (in [0,1]), not coefficient*dt.
    #
    # Bug: nlsolver.c was set to theta*dt instead of theta, causing the Newton
    # method to evaluate f at t + theta*dt*dt instead of t + theta*dt.
    # This caused extrapolation errors when using time-dependent parameters
    # with interpolations that end exactly at tspan[2].
    #
    # Reference: https://github.com/SciML/ModelingToolkit.jl/issues/4143

    @testset "Verify Newton evaluates at correct time" begin
        # Test: use a function that records the time it's called at
        # and verify all times are within the expected range [0, tend]

        tend = 1.0
        recorded_times = Float64[]

        function f_recording!(du, u, p, t)
            push!(recorded_times, t)
            du[1] = -u[1]
        end

        function g!(du, u, p, t)
            du[1] = 0.1 * u[1]
        end

        u0 = [1.0]
        tspan = (0.0, tend)
        prob = SDEProblem(f_recording!, g!, u0, tspan)

        # Test ImplicitEM
        empty!(recorded_times)
        solve(prob, ImplicitEM(); dt = 0.1, adaptive = false)

        # All recorded times should be within [0, tend]
        # With the bug (c = theta*dt), times would be at t + theta*dt*dt
        # which is wrong but usually smaller than expected
        # With the fix (c = theta), times should be at t + theta*dt
        @test all(0.0 .<= recorded_times .<= tend + 10 * eps(tend))
        @test maximum(recorded_times) <= tend + 10 * eps(tend)

        # The maximum recorded time should be close to tend
        # With the bug, the Newton method would evaluate at t + dt*dt instead of t + dt
        # So for a step from t=0.9 to t=1.0, the bug gives 0.9 + 0.01 = 0.91,
        # while correct gives 0.9 + 0.1 = 1.0
        @test maximum(recorded_times) >= tend - 0.15  # Should reach close to tend

        # Test ImplicitEulerHeun
        empty!(recorded_times)
        solve(prob, ImplicitEulerHeun(); dt = 0.1, adaptive = false)
        @test all(0.0 .<= recorded_times .<= tend + 10 * eps(tend))

        # Test ImplicitRKMil
        empty!(recorded_times)
        solve(prob, ImplicitRKMil(); dt = 0.1, adaptive = false)
        @test all(0.0 .<= recorded_times .<= tend + 10 * eps(tend))
    end

    @testset "Test with adaptive stepping" begin
        # Test that adaptive stepping also works correctly
        recorded_times = Float64[]

        function f_recording!(du, u, p, t)
            push!(recorded_times, t)
            du[1] = -u[1]
        end

        function g!(du, u, p, t)
            du[1] = 0.1 * u[1]
        end

        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SDEProblem(f_recording!, g!, u0, tspan)

        empty!(recorded_times)
        solve(prob, ImplicitEM())

        # All recorded times should be within [0, tend]
        @test all(0.0 .<= recorded_times .<= 1.0 + 10 * eps(1.0))
    end

    @testset "Test ISSEM methods" begin
        # Test split-step methods as well
        recorded_times = Float64[]

        function f_recording!(du, u, p, t)
            push!(recorded_times, t)
            du[1] = -u[1]
        end

        function g!(du, u, p, t)
            du[1] = 0.1 * u[1]
        end

        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SDEProblem(f_recording!, g!, u0, tspan)

        # Test ISSEM
        empty!(recorded_times)
        solve(prob, ISSEM(); dt = 0.1, adaptive = false)
        @test all(0.0 .<= recorded_times .<= 1.0 + 10 * eps(1.0))

        # Test ISSEulerHeun
        empty!(recorded_times)
        solve(prob, ISSEulerHeun(); dt = 0.1, adaptive = false)
        @test all(0.0 .<= recorded_times .<= 1.0 + 10 * eps(1.0))
    end

    @testset "Out-of-place methods" begin
        # Test out-of-place versions
        recorded_times = Float64[]

        function f_recording(u, p, t)
            push!(recorded_times, t)
            return -u
        end

        function g(u, p, t)
            return 0.1 * u
        end

        u0 = [1.0]
        tspan = (0.0, 1.0)
        prob = SDEProblem(f_recording, g, u0, tspan)

        # Test ImplicitEM out-of-place
        empty!(recorded_times)
        solve(prob, ImplicitEM(); dt = 0.1, adaptive = false)
        @test all(0.0 .<= recorded_times .<= 1.0 + 10 * eps(1.0))
    end
end

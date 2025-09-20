using OrdinaryDiffEq, Test, Random, StaticArrays, DiffEqCallbacks
import ODEProblemLibrary: prob_ode_linear
Random.seed!(100)
@testset "Tstops Tests on the Interval [0, 1]" begin
    prob = prob_ode_linear

    sol = solve(prob, Tsit5(), dt = 1 // 2^(6), tstops = [1 / 2])
    @test 1 // 2 ∈ sol.t

    sol = solve(prob, RK4(), dt = 1 // 3, tstops = [1 / 2], adaptive = false)
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    sol = solve(
        prob, RK4(), dt = 1 // 3, tstops = [1 / 2],
        d_discontinuities = [-1 / 2, 1 / 2, 3 / 2], adaptive = false
    )
    @test sol.t == [0, 1 / 3, 1 / 2, 1 / 3 + 1 / 2, 1]

    # TODO
    integrator = init(
        prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4],
        adaptive = false
    )

    sol = solve(prob, RK4(), tstops = [1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4], adaptive = false)
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(
        prob, RK4(), tstops = [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1],
        adaptive = false
    )
    @test sol.t == [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]

    sol = solve(prob, RK4(), tstops = 0:(1 // 16):1, adaptive = false)
    @test sol.t == collect(0:(1 // 16):1)

    sol = solve(prob, RK4(), tstops = range(0, stop = 1, length = 100), adaptive = false)
    @test sol.t == collect(range(0, stop = 1, length = 100))
end

@testset "Integrator Tstops Tests on the Interval $(["[-1, 0]", "[0, 1]"][i])" for (i, tdir) in enumerate(
        [
            -1.0;
            1.0
        ]
    )
    prob2 = remake(prob_ode_linear, tspan = (0.0, tdir * 1.0))
    integrator = init(prob2, Tsit5())
    tstops = tdir .* [0, 1 / 5, 1 / 4, 1 / 3, 1 / 2, 3 / 4, 1]
    for tstop in tstops
        add_tstop!(integrator, tstop)
    end
    @test_throws ErrorException add_tstop!(integrator, -0.1 * tdir)
    solve!(integrator)
    for tstop in tstops
        @test tstop ∈ integrator.sol.t
    end
end

@testset "Tstops Eps" begin
    function de(du, u, p, t) # specific DE does not impact the issue
        a, b = p
        du[1] = a * u[1]
        du[2] = b * u[2]
    end

    saveat = [0.0, 0.0094777, 1.5574]
    tstop = 0.010823

    affect!(integrator) = integrator.u[1] += 1.0
    condition(u, t, integrator) = t == tstop
    callback = DiscreteCallback(condition, affect!)

    prob = ODEProblem(de, zeros(2), (-1, 3.0), rand(2))
    sol = solve(prob, Tsit5(), saveat = saveat, tstops = tstop, callback = callback)
    @test sol.t[end] == 1.5574
end

@testset "Tstops Type Conversion" begin
    called = Ref(false)
    tval = rand()
    ff(du, u, p, t) = du .= 0
    cb = DiscreteCallback(
        (u, t, integrator) -> t == Float32(tval),
        integrator -> (called[] = true)
    )
    prob = ODEProblem(ff, [0.0], (0.0f0, 1.0f0))
    sol = solve(prob, Tsit5(), tstops = [tval], callback = cb)
end

@testset "Late binding tstops" begin
    function rhs(u, p, t)
        u * p + t
    end
    prob = ODEProblem(rhs, 1.0, (0.0, 1.0), 0.1; tstops = (p, tspan) -> tspan[1]:p:tspan[2])
    sol = solve(prob, Tsit5())
    @test 0.0:0.1:1.0 ⊆ sol.t
    prob2 = remake(prob; p = 0.07)
    sol2 = solve(prob2, Tsit5())
    @test 0.0:0.07:1.0 ⊆ sol2.t
end

@testset "Tstop Overshoot and Dense Time Event Tests" begin
    # Tests for issue #2752: tstop overshoot errors with StaticArrays
    
    @testset "StaticArrays vs Arrays with extreme precision" begin
        # Test the specific case that was failing: extreme precision + StaticArrays
        function precise_dynamics(u, p, t)
            x = @view u[1:2]
            v = @view u[3:4]
            
            # Electromagnetic-like dynamics
            dv = -0.01 * x + 1e-6 * sin(100*t) * SVector{2}(1, 1)
            
            return SVector{4}(v[1], v[2], dv[1], dv[2])
        end
        
        function precise_dynamics_array!(du, u, p, t)
            x = @view u[1:2]
            v = @view u[3:4]
            
            dv = -0.01 * x + 1e-6 * sin(100*t) * [1, 1]
            du[1] = v[1]
            du[2] = v[2]  
            du[3] = dv[1]
            du[4] = dv[2]
        end
        
        # Initial conditions
        u0_static = SVector{4}(1.0, -0.5, 0.01, 0.01)
        u0_array = [1.0, -0.5, 0.01, 0.01]
        tspan = (0.0, 2.0)
        tstops = [0.5, 1.0, 1.5]
        
        # Test with extreme tolerances that originally caused issues
        prob_static = ODEProblem(precise_dynamics, u0_static, tspan)
        sol_static = solve(prob_static, Vern9(); reltol=1e-12, abstol=1e-15, 
                          tstops=tstops)
        @test SciMLBase.successful_retcode(sol_static)
        for tstop in tstops
            @test tstop ∈ sol_static.t
        end
        
        prob_array = ODEProblem(precise_dynamics_array!, u0_array, tspan)
        sol_array = solve(prob_array, Vern9(); reltol=1e-12, abstol=1e-15, 
                         tstops=tstops)
        @test SciMLBase.successful_retcode(sol_static)
        for tstop in tstops
            @test tstop ∈ sol_array.t
        end
        
        # Solutions should be very close despite different array types
        @test isapprox(sol_static(2.0), sol_array(2.0), rtol=1e-10)
    end
    
    @testset "Duplicate tstops handling" begin
        function simple_ode(u, p, t)
            SA[0.1 * u[1]]
        end
        
        u0 = SVector{1}(1.0)
        tspan = (0.0, 2.0)
        
        # Test multiple identical tstops - should all be processed
        duplicate_tstops = [0.5, 0.5, 0.5, 1.0, 1.0]
        
        prob = ODEProblem(simple_ode, u0, tspan)
        sol = solve(prob, Vern9(); tstops=duplicate_tstops)
        
        @test SciMLBase.successful_retcode(sol)
        
        # Count how many times each tstop appears in solution
        count_05 = count(t -> abs(t - 0.5) < 1e-12, sol.t)
        count_10 = count(t -> abs(t - 1.0) < 1e-12, sol.t)
        
        # Should handle all duplicate tstops (though may not save all due to deduplication)
        @test count_05 >= 1  # At least one 0.5
        @test count_10 >= 1  # At least one 1.0
        
        # Test with StaticArrays too
        prob_static = ODEProblem(simple_ode, u0, tspan)
        sol_static = solve(prob_static, Vern9(); tstops=duplicate_tstops)
        @test SciMLBase.successful_retcode(sol_static)
    end
    
    @testset "PresetTimeCallback with identical times" begin
        # Test PresetTimeCallback scenarios where callbacks are set at same times as tstops
        
        event_times = Float64[]
        callback_times = Float64[]
        
        function affect_preset!(integrator)
            push!(callback_times, integrator.t)
            integrator.u += 0.1* integrator.u  # Small modification
        end
        
        function simple_growth(u, p, t)
            SA[0.1 * u[1]]
        end
        
        u0 = SA[1.0]
        tspan = (0.0, 3.0)
        
        # Define times where both tstops and callbacks should trigger
        critical_times = [0.5, 1.0, 1.5, 2.0, 2.5]
        
        # Create PresetTimeCallback at the same times as tstops
        preset_cb = PresetTimeCallback(critical_times, affect_preset!)
        
        prob = ODEProblem(simple_growth, u0, tspan)
        sol = solve(prob, Vern9(); tstops=critical_times, callback=preset_cb, 
                   reltol=1e-10, abstol=1e-12)
        
        @test SciMLBase.successful_retcode(sol)
        
        # Verify all tstops were hit
        for time in critical_times
            @test any(abs.(sol.t .- time) .< 1e-10)
        end
        
        # Verify all callbacks were triggered  
        @test length(callback_times) == length(critical_times)
        for time in critical_times
            @test any(abs.(callback_times .- time) .< 1e-10)
        end
        
        # Test the same with regular arrays
        u0_array = [1.0]
        callback_times_array = Float64[]
        
        function affect_preset_array!(integrator)
            push!(callback_times_array, integrator.t)
            integrator.u[1] += 0.1
        end
        
        function simple_growth_array!(du, u, p, t)
            du[1] = 0.1 * u[1]
        end
        
        preset_cb_array = PresetTimeCallback(critical_times, affect_preset_array!)
        
        prob_array = ODEProblem(simple_growth_array!, u0_array, tspan)
        sol_array = solve(prob_array, Vern9(); tstops=critical_times, callback=preset_cb_array,
                         reltol=1e-10, abstol=1e-12)
        
        @test SciMLBase.successful_retcode(sol_array)
        @test length(callback_times_array) == length(critical_times)
        
        # Both should have triggered all events
        @test length(callback_times) == length(callback_times_array) == length(critical_times)
    end
    
    @testset "Tiny tstop step handling" begin
        # Test cases where tstop is very close to current time (dt < eps(t))
        function test_ode(u, p, t)
            SA[0.01 * u[1]]
        end
        
        u0 = SVector{1}(1.0)
        tspan = (0.0, 1.0)
        
        # Create tstop very close to start time (would cause tiny dt)
        tiny_tstops = [1e-15, 1e-14, 1e-13]
        
        for tiny_tstop in tiny_tstops
            prob = ODEProblem(test_ode, u0, tspan)
            sol = solve(prob, Vern9(); tstops=[tiny_tstop])
            @test SciMLBase.successful_retcode(sol)
            @test any(abs.(sol.t .- tiny_tstop) .< 1e-14)  # Should handle tiny tstop correctly
        end
    end
    
    @testset "Multiple close tstops with StaticArrays" begin
        # Test with multiple tstops that are very close together - stress test the flag logic
        function oscillator(u, p, t)
            SVector{2}(u[2], -u[1])  # Simple harmonic oscillator
        end
        
        u0 = SVector{2}(1.0, 0.0)
        tspan = (0.0, 4.0)
        
        # Multiple tstops close together (within floating-point precision range)
        close_tstops = [1.0, 1.0 + 1e-14, 1.0 + 2e-14, 1.0 + 5e-14, 
                       2.0, 2.0 + 1e-15, 2.0 + 1e-14,
                       3.0, 3.0 + 1e-13]
        
        prob = ODEProblem(oscillator, u0, tspan)
        sol = solve(prob, Vern9(); tstops=close_tstops, reltol=1e-12, abstol=1e-15)
        
        @test SciMLBase.successful_retcode(sol)
        
        # Should handle all close tstops without error
        # (Some might be deduplicated, but no errors should occur)
        unique_times = [1.0, 2.0, 3.0]
        for time in unique_times
            @test any(abs.(sol.t .- time) .< 1e-10)  # At least hit the main times
        end
    end
    
    @testset "Backward integration with tstop flags" begin
        # Test that the fix works for backward time integration
        function decay_ode(u, p, t)
            SA[-0.1 * u[1]]
        end
        
        u0 = SVector{1}(1.0)
        tspan = (2.0, 0.0)  # Backward integration
        tstops = [1.5, 1.0, 0.5]
        
        prob = ODEProblem(decay_ode, u0, tspan)
        sol = solve(prob, Vern9(); tstops=tstops, reltol=1e-12, abstol=1e-15)
        
        @test SciMLBase.successful_retcode(sol)
        for tstop in tstops
            @test tstop ∈ sol.t
        end
    end
    
    @testset "Continuous callbacks during tstop steps" begin
        # Test that continuous callbacks work properly with tstop flag mechanism
        
        crossing_times = Float64[]
        
        function affect_continuous!(integrator)
            push!(crossing_times, integrator.t)
        end
        
        function condition_continuous(u, t, integrator)
            u[1] - 0.5  # Crosses when u[1] = 0.5
        end
        
        function exponential_growth(u, p, t)
            [0.2 * u[1]]  # Exponential growth
        end
        
        u0 = [0.1]  # Start below 0.5
        tspan = (0.0, 10.0)
        tstops = [2.0, 4.0, 6.0, 8.0]  # Regular tstops
        
        continuous_cb = ContinuousCallback(condition_continuous, affect_continuous!)
        
        prob = ODEProblem(exponential_growth, u0, tspan)
        sol = solve(prob, Vern9(); tstops=tstops, callback=continuous_cb,
                   reltol=1e-10, abstol=1e-12)
        
        @test SciMLBase.successful_retcode(sol)
        
        # Should hit all tstops
        for tstop in tstops
            @test tstop ∈ sol.t
        end
        
        # Should also detect continuous callback crossings
        @test length(crossing_times) > 0  # At least one crossing detected
        
        # Verify crossings are at correct value
        for crossing_time in crossing_times
            u_at_crossing = sol(crossing_time)
            @test abs(u_at_crossing[1] - 0.5) < 1e-8  # Should be very close to 0.5
        end
    end
    
end

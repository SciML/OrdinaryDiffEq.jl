using OrdinaryDiffEqVerner, StaticArrays, Test
using OrdinaryDiffEqCore: handle_tstop_step!

# Test cases for the tstop robustness fix with next_step_tstop flag
@testset "Tstop Robustness Tests" begin
    
    @testset "Basic tstop flag functionality" begin
        # Simple ODE problem to test flag behavior
        function simple_ode!(du, u, p, t)
            du[1] = -u[1]
            du[2] = u[1] - u[2]
        end
        
        function simple_ode(u, p, t)
            [-u[1], u[1] - u[2]]
        end
        
        u0_array = [1.0, 0.0]
        u0_static = SVector{2}(1.0, 0.0)
        tspan = (0.0, 1.0)
        
        # Test with regular arrays
        prob_array = ODEProblem(simple_ode!, u0_array, tspan)
        sol_array = solve(prob_array, Vern9(); reltol=1e-12, abstol=1e-12, 
                         tstops=[0.5], save_everystep=false)
        @test sol_array.retcode == :Success
        @test 0.5 in sol_array.t  # Should have saved at tstop
        
        # Test with StaticArrays - should work now without tstop error
        prob_static = ODEProblem(simple_ode, u0_static, tspan)
        sol_static = solve(prob_static, Vern9(); reltol=1e-12, abstol=1e-12, 
                          tstops=[0.5], save_everystep=false)
        @test sol_static.retcode == :Success
        @test 0.5 in sol_static.t  # Should have saved at tstop
        
        # Solutions should be very close despite different array types
        @test isapprox(sol_array(1.0), sol_static(1.0), rtol=1e-10)
    end
    
    @testset "Tiny tstop step handling" begin
        # Test case where tstop is very close to current time
        function test_ode(u, p, t)
            [u[1]]  # Simple growth
        end
        
        u0 = SVector{1}(1.0)
        tspan = (0.0, 1.0)
        
        # Create tstop very close to start time (would cause tiny dt)
        tiny_tstop = 1e-15
        
        prob = ODEProblem(test_ode, u0, tspan)
        sol = solve(prob, Vern9(); tstops=[tiny_tstop], save_everystep=false)
        
        @test sol.retcode == :Success
        @test tiny_tstop in sol.t  # Should handle tiny tstop correctly
    end
    
    @testset "Multiple close tstops" begin
        # Test with multiple tstops that are very close together
        function growth_ode(u, p, t)
            [0.1 * u[1]]
        end
        
        u0 = SVector{1}(1.0)
        tspan = (0.0, 2.0)
        
        # Multiple tstops close together
        close_tstops = [0.5, 0.5 + 1e-14, 0.5 + 2e-14, 1.0]
        
        prob = ODEProblem(growth_ode, u0, tspan)
        sol = solve(prob, Vern9(); tstops=close_tstops, reltol=1e-12, abstol=1e-12)
        
        @test sol.retcode == :Success
        # All tstops should be handled correctly
        for tstop in close_tstops
            @test any(abs.(sol.t .- tstop) .< 1e-12)  # Should have hit each tstop
        end
    end
    
    @testset "Extreme precision with StaticArrays" begin
        # Test the specific case that was failing: extreme precision + StaticArrays
        function precise_dynamics(u, p, t)
            # Simplified electromagnetic-like dynamics
            x = @view u[1:2]
            v = @view u[3:4]
            
            # Simple force model
            dv = -0.01 * x + 1e-6 * sin(1000*t) * [1, 1]
            
            return SVector{4}(v[1], v[2], dv[1], dv[2])
        end
        
        # Initial conditions similar to the original issue
        u0 = SVector{4}(1.0, -0.5, 0.01, 0.01)
        tspan = (-1.0, 1.0)
        
        # Test with extreme tolerances that originally caused issues
        prob = ODEProblem(precise_dynamics, u0, tspan)
        sol = solve(prob, Vern9(); reltol=1e-12, abstol=1e-15, 
                   tstops=[0.0], save_everystep=false, maxiters=10^6)
        
        @test sol.retcode == :Success
        @test 0.0 in sol.t
    end
    
    @testset "Flag state management" begin
        # Test that flags are properly set and reset
        function flag_test_ode(u, p, t)
            [u[1]]
        end
        
        u0 = SVector{1}(1.0)
        prob = ODEProblem(flag_test_ode, u0, (0.0, 2.0))
        
        # Create integrator manually to inspect flag states
        integrator = init(prob, Vern9(); tstops=[1.0])
        
        # Initially, flag should be false
        @test integrator.next_step_tstop == false
        
        # Step until we approach the tstop
        while integrator.t < 0.9
            step!(integrator)
            # Flag should still be false when not near tstop
            @test integrator.next_step_tstop == false
        end
        
        # Take steps near tstop - flag should get set
        while integrator.t < 1.0
            step!(integrator)
            # When dt is reduced for tstop, flag should be set
            if integrator.next_step_tstop
                @test integrator.tstop_target ≈ 1.0
                break
            end
        end
        
        # After hitting tstop, flag should be reset
        step!(integrator)
        @test integrator.next_step_tstop == false
        
        finalize!(integrator)
    end
    
    @testset "Backward time integration" begin
        # Test that the fix works for backward time integration too
        function backward_ode(u, p, t)
            [-u[1]]  # Decay
        end
        
        u0 = SVector{1}(1.0)
        tspan = (1.0, 0.0)  # Backward integration
        
        prob = ODEProblem(backward_ode, u0, tspan)
        sol = solve(prob, Vern9(); tstops=[0.5], reltol=1e-12, abstol=1e-12)
        
        @test sol.retcode == :Success
        @test 0.5 in sol.t
    end
    
end
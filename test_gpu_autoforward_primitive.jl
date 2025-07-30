"""
Test for GPU-safe AutoForwardFromPrimitive wrapper functionality.

This test verifies that OrdinaryDiffEq automatically wraps AutoForwardDiff 
with AutoForwardFromPrimitive when dealing with GPU arrays to avoid scalar indexing issues.
"""

using OrdinaryDiffEq, JLArrays, LinearAlgebra, ADTypes, Test

function f(du, u, p, t)
    A = jl(-ones(3, 3))
    return mul!(du, A, u)
end

@testset "GPU AutoForwardFromPrimitive Tests" begin
    u0 = jl([1.0; 0.0; 0.0])
    tspan = (0.0f0, 10.0f0)
    prob = ODEProblem(f, u0, tspan)

    @testset "Rosenbrock23 with AutoForwardDiff on GPU arrays" begin
        # This should now work without scalar indexing errors due to automatic wrapping
        sol = solve(prob, Rosenbrock23(; autodiff=AutoForwardDiff()))
        @test sol.retcode == ReturnCode.Success
        @test length(sol.u) > 1
        
        # Verify the solution has reasonable values
        @test all(isfinite.(sol[end]))
    end

    @testset "Compare GPU vs CPU solutions" begin
        # GPU solution 
        gpu_sol = solve(prob, Rosenbrock23(; autodiff=AutoForwardDiff()))
        
        # CPU equivalent
        u0_cpu = [1.0, 0.0, 0.0]
        function f_cpu(du, u, p, t)
            A = -ones(3, 3)
            return mul!(du, A, u)
        end
        prob_cpu = ODEProblem(f_cpu, u0_cpu, tspan)
        cpu_sol = solve(prob_cpu, Rosenbrock23(; autodiff=AutoForwardDiff()))
        
        # Solutions should be very similar
        @test isapprox(Array(gpu_sol[end]), cpu_sol[end], rtol=1e-10)
    end

    @testset "Different Rosenbrock methods" begin
        # Test other Rosenbrock methods that use autodiff
        for alg in [Rosenbrock23(), Rodas4(), Rodas5P()]
            sol = solve(prob, alg)
            @test sol.retcode == ReturnCode.Success
        end
    end
end

println("✅ All GPU AutoForwardFromPrimitive tests passed!")
println("✅ OrdinaryDiffEq now automatically handles GPU arrays with AutoForwardDiff!")
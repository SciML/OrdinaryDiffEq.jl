using OrdinaryDiffEq, JLArrays, LinearAlgebra, Test, ADTypes

@testset "GPU AutoDiff with JLArrays" begin
    function f(du, u, p, t)
        A = jl(-ones(3, 3))
        return mul!(du, A, u)
    end

    u0 = jl([1.0; 0.0; 0.0])
    tspan = (0.0f0, 100.0f0)
    prob = ODEProblem(f, u0, tspan)

    # This should fail with scalar indexing error when using AutoForwardDiff (default)
    @test_throws Exception solve(prob, Rosenbrock23())
end
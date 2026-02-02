using OrdinaryDiffEq, JLArrays, LinearAlgebra, Test, ADTypes

@testset "GPU AutoDiff with JLArrays" begin
    function f(u, p, t)
        A = jl(-ones(3, 3))
        return A * u
    end
    function f!(du, u, p, t)
        A = jl(-ones(3, 3))
        return mul!(du, A, u)
    end

    u0 = jl([1.0; 0.0; 0.0])
    tspan = (0.0f0, 100.0f0)
    prob = ODEProblem{false}(f, u0, tspan)
    solve(prob, TRBDF2())
    solve(prob, Rodas5P())

    prob2 = ODEProblem(f!, u0, tspan)
    solve(prob, TRBDF2())
    solve(prob, Rodas5P())
end

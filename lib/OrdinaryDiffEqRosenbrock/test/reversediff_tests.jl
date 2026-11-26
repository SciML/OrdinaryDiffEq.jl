using OrdinaryDiffEqRosenbrock, ReverseDiff, Test
using ADTypes: AutoFiniteDiff
using DiffEqBase: SensitivityADPassThrough

function reverse_rhs!(du, u, p, t)
    du[1] = -p[1] * u[1]
    du[2] = p[1] * u[1] - u[2]
    return nothing
end

function reverse_loss(p)
    prob = ODEProblem(reverse_rhs!, eltype(p).([1.0, 0.0]), (0.0, 1.0), p)
    sol = solve(
        prob, Rodas5(autodiff = AutoFiniteDiff());
        sensealg = SensitivityADPassThrough(), abstol = 1.0e-6, reltol = 1.0e-6
    )
    return sol.u[end][1]
end

@testset "ReverseDiff through Rodas5" begin
    @test_throws ArgumentError("ReverseDiff through an in-place Rosenbrock solve with AutoFiniteDiff and the default linear solver is unsupported. Use ForwardDiff or a SciMLSensitivity sensitivity algorithm instead of SensitivityADPassThrough().") ReverseDiff.gradient(reverse_loss, [2.0])
end

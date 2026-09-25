using OrdinaryDiffEqCore
using Test

struct VectorSolverCache
    nlsolver::Vector{Int}
end

struct FakeIntegrator
    cache::VectorSolverCache
end

@testset "instability diagnostics with a vector of nonlinear solvers" begin
    integrator = FakeIntegrator(VectorSolverCache([1, 2]))
    @test OrdinaryDiffEqCore.instability_jacobian(integrator) === nothing
end

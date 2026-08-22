using OrdinaryDiffEqCore
using Test

# `instability_jacobian` reaches for `cache.nlsolver.cache.J` to report the
# Jacobian behind a blow-up. Caches that run several nonlinear solves in
# parallel store a vector of solvers in that field, and reaching through it
# threw instead of reporting a diagnostic.
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

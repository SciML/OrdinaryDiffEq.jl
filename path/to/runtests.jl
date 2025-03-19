using Test
using OrdinaryDiffEqTaylorSeries: build_index_reduced_system

# A small dummy DAE function:
function dummy_dae(du, x, p, t)
    # Not doing actual computations - just a placeholder
    return du .+ x
end

# A mock-up compute_system_jacobian for testing:
function compute_system_jacobian(f, taylor_u, p, t0, order)
    # Return a trivial 2x2 identity matrix so we can confirm the function gets called.
    return [1.0 0.0; 0.0 1.0]
end

@testset "DAETS build_index_reduced_system test" begin
    # Let's pick a small initial condition
    u0 = [1.0, 2.0]
    p = nothing
    t0 = 0.0

    # We call our function under test:
    order, jac = build_index_reduced_system(dummy_dae, u0, p, t0; max_order=5)

    @test order == 0  # Because our identity matrix is already nonsingular
    @test jac == [1.0 0.0; 0.0 1.0]
end 
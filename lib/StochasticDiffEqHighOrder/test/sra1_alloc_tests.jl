using StochasticDiffEqHighOrder
using OrdinaryDiffEqCore
using SciMLBase
using Test

# Discard return values so @allocated does not measure ReturnCode boxing (Julia 1.10).
@inline function perform_step_void!(integrator, cache)
    OrdinaryDiffEqCore.perform_step!(integrator, cache)
    return nothing
end

@testset "SRA1 iip k₁ *= dt is in-place (no per-step alloc)" begin
    # Elementwise rhs (same style as StochasticDiffEq alloc_tests) so AutoSpecialize
    # does not inject unrelated broadcast temporaries into the measurement.
    f!(du, u, p, t) = (du[1] = -u[1]; du[2] = -u[2]; du[3] = -u[3]; nothing)
    g!(du, u, p, t) = (du[1] = 0.1; du[2] = 0.2; du[3] = 0.3; nothing)
    u0 = [1.0, 1.0, 1.0]
    prob = SDEProblem(f!, g!, u0, (0.0, 10.0))

    # perform_step! isolates the k₁ *= dt rebinding (96 B/step on master for n=4).
    integ = init(prob, SRA1(); adaptive = false, dt = 0.01, save_on = false)
    cache = integ.cache
    for _ in 1:10
        perform_step_void!(integ, cache)
        integ.t = mod(integ.t, 9.0)
    end
    @test @allocated(perform_step_void!(integ, cache)) == 0
end

@testset "SRA1 iip/oop agree after in-place k₁ scale" begin
    f!(du, u, p, t) = (du[1] = -u[1]; du[2] = -u[2]; du[3] = -u[3]; du[4] = -u[4]; nothing)
    f(u, p, t) = -u
    g!(du, u, p, t) = (du[1] = 0.1; du[2] = 0.1; du[3] = 0.1; du[4] = 0.1; nothing)
    g(u, p, t) = fill(0.1, length(u))
    u0 = [1.0, 0.5, 0.25, 0.125]
    tspan = (0.0, 1.0)
    prob_iip = SDEProblem(f!, g!, u0, tspan)
    prob_oop = SDEProblem(f, g, u0, tspan)
    sol_iip = solve(prob_iip, SRA1(); seed = 7, adaptive = false, dt = 1 / 64)
    sol_oop = solve(prob_oop, SRA1(); seed = 7, adaptive = false, dt = 1 / 64)
    @test sol_iip.u[end] ≈ sol_oop.u[end] rtol = 2.0e-16 atol = 2.0e-16
end

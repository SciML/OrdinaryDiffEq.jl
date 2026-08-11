using OrdinaryDiffEqRosenbrock, Test

@testset "Rodas5P initial-time interpolation" begin
    prob = ODEProblem((u, p, t) -> -u, [1.0], (0.0, 1.0))
    integrator = init(prob, Rodas5P())
    t0 = first(integrator.sol.t)

    @test integrator.sol(t0) == prob.u0
    @test only(integrator.sol([t0]).u) == prob.u0
    @test integrator.sol(t0; idxs = 1) == only(prob.u0)
    @test only(integrator.sol([t0]; idxs = 1).u) == only(prob.u0)

    out = similar(prob.u0)
    integrator.sol(out, t0)
    @test out == prob.u0

    outs = [similar(prob.u0)]
    integrator.sol(outs, [t0])
    @test only(outs) == prob.u0

    scalar_outs = zeros(1)
    integrator.sol(scalar_outs, [t0]; idxs = 1)
    @test only(scalar_outs) == only(prob.u0)
end

using OrdinaryDiffEqBDF, OrdinaryDiffEqSymplecticRK, OrdinaryDiffEqFunctionMap, Test

@testset "Default dense output" begin
    res(du, u, p, t) = [du[1] + u[1], u[1] + u[2] - 1]
    dprob = DAEProblem(res, [-1.0, 0.0], [1.0, 0.0], (0.0, 1.0); differential_vars = [true, false])
    @test !solve(dprob, DImplicitEuler()).dense
    @test !solve(dprob, DABDF2()).dense
    @test solve(dprob, DFBDF()).dense
    @test solve(dprob, DNordsieckBDF()).dense

    sprob = SecondOrderODEProblem((v, u, p, t) -> -u, [0.0], [1.0], (0.0, 1.0))
    @test !solve(sprob, VerletLeapfrog(); dt = 0.1).dense
    @test !solve(sprob, LeapfrogDriftKickDrift(); dt = 0.1).dense
    @test solve(sprob, PseudoVerletLeapfrog(); dt = 0.1).dense

    @test !solve(DiscreteProblem((u, p, t) -> 1.1u, [1.0], (0.0, 10.0)), FunctionMap()).dense
end

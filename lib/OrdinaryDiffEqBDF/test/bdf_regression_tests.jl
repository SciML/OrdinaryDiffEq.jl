using OrdinaryDiffEqBDF, Test

foop = (u, p, t) -> u
proboop = ODEProblem(foop, ones(2), (0.0, 1000.0))

fiip = (du, u, p, t) -> du .= u
probiip = ODEProblem(fiip, ones(2), (0.0, 1000.0))

@testset "FBDF reinit" begin
    for prob in [proboop, probiip]
        integ = init(prob, FBDF(), verbose = false) #suppress warning to clean up CI
        solve!(integ)
        @test integ.sol.retcode != ReturnCode.Success
        @test integ.sol.t[end] >= 700
        reinit!(integ, prob.u0)
        solve!(integ)
        @test integ.sol.retcode != ReturnCode.Success
        @test integ.sol.t[end] >= 700
    end
end

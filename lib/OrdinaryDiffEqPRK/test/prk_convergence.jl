# This definitely needs cleaning
using OrdinaryDiffEqPRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

dts2 = 1 .// 2 .^ (7:-1:3)
testTol = 0.2

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]
    sim3 = test_convergence(dts2, prob, KuttaPRK2p5(threading = true))
    @test sim3.𝒪est[:l∞] ≈ 5 atol = testTol

    sim3 = test_convergence(dts2, prob, KuttaPRK2p5(threading = false))
    @test sim3.𝒪est[:l∞] ≈ 5 atol = testTol
end

@testset "Threaded stages inside another threaded loop" begin
    prob_iip = ODEProblem((du, u, p, t) -> (du .= -u), [1.0, 2.0], (0.0, 1.0))
    prob_oop = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    for prob in (prob_iip, prob_oop)
        serial = solve(prob, KuttaPRK2p5(threading = false); dt = 0.1)
        sols = Vector{Any}(undef, 2)
        Threads.@threads for i in 1:2
            sols[i] = solve(prob, KuttaPRK2p5(threading = true); dt = 0.1)
        end
        @test all(sol -> sol.u == serial.u, sols)
    end
end

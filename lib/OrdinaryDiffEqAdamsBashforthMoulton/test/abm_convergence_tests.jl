# This definitely needs cleaning
using OrdinaryDiffEqAdamsBashforthMoulton, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

## Convergence Testing
dts1 = 1 .// 2 .^ (9:-1:5)
dts = 1 .// 2 .^ (8:-1:4)
testTol = 0.2

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in
    1:2

    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    sim5 = test_convergence(dts, prob, AB3())
    @test sim5.𝒪est[:l2] ≈ 3 atol = testTol
    sim6 = test_convergence(dts, prob, ABM32())
    @test sim6.𝒪est[:l2] ≈ 3 atol = testTol
    sim7 = test_convergence(dts, prob, AB4())
    @test sim7.𝒪est[:l2] ≈ 4 atol = testTol
    sim8 = test_convergence(dts1, prob, ABM43())  #using dts1 due to floating point error in convergence test
    @test sim8.𝒪est[:l2] ≈ 4 atol = testTol
    sim9 = test_convergence(dts, prob, AB5())
    @test sim9.𝒪est[:l2] ≈ 5 atol = testTol
    sim10 = test_convergence(dts, prob, ABM54())
    @test sim10.𝒪est[:l2] ≈ 5 atol = testTol
    sim101 = test_convergence(dts, prob, VCAB3())
    @test sim101.𝒪est[:l2] ≈ 3 atol = testTol
    sim102 = test_convergence(dts, prob, VCAB4())
    @test sim102.𝒪est[:l2] ≈ 4 atol = testTol
    sim103 = test_convergence(dts, prob, VCAB5())
    @test sim103.𝒪est[:l2] ≈ 5 atol = testTol
    sim104 = test_convergence(dts, prob, VCABM3())
    @test sim104.𝒪est[:l2] ≈ 3 atol = testTol
    sim105 = test_convergence(dts, prob, VCABM4())
    @test sim105.𝒪est[:l2] ≈ 4 atol = testTol
    sim106 = test_convergence(dts, prob, VCABM5())
    @test sim106.𝒪est[:l2] ≈ 5 atol = testTol
end

@testset "ABM steady-step function evaluations" begin
    for (alg, startup) in ((ABM32(), 2), (ABM43(), 3), (ABM54(), 4))
        for inplace in (false, true)
            calls = Ref(0)
            prob = if inplace
                ODEProblem(
                    (du, u, p, t) -> (calls[] += 1; du[1] = -u[1]),
                    [1.0], (0.0, 1.0)
                )
            else
                ODEProblem((u, p, t) -> (calls[] += 1; -u), 1.0, (0.0, 1.0))
            end
            integrator = init(prob, alg; dt = 1 / 32, adaptive = false)
            for _ in 1:startup
                step!(integrator)
            end
            for _ in 1:4
                nf = integrator.stats.nf
                old_calls = calls[]
                step!(integrator)
                @test integrator.stats.nf - nf == 2
                @test calls[] - old_calls == 2
            end
        end
    end
end

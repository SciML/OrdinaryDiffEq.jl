using OrdinaryDiffEqAdamsBashforthMoulton, ODEProblemLibrary, DiffEqDevTools
using FastBroadcast: Threaded
using Polyester
using Test, Random
Random.seed!(100)

dts = 1 .// 2 .^ (8:-1:4)
testTol = 0.2

# Threaded broadcasting via Polyester only supports vectors, not matrices
@testset "Explicit Solver Convergence Tests (out-of-place) - threaded" begin
    prob = ODEProblemLibrary.prob_ode_linear

    sim5 = test_convergence(dts, prob, AB3(thread = Threaded()))
    @test sim5.𝒪est[:l2] ≈ 3 atol = testTol
    sim7 = test_convergence(dts, prob, AB4(thread = Threaded()))
    @test sim7.𝒪est[:l2] ≈ 4 atol = testTol
    sim9 = test_convergence(dts, prob, AB5(thread = Threaded()))
    @test sim9.𝒪est[:l2] ≈ 5 atol = testTol
    sim101 = test_convergence(dts, prob, VCAB3(thread = Threaded()))
    @test sim101.𝒪est[:l2] ≈ 3 atol = testTol
    sim103 = test_convergence(dts, prob, VCAB5(thread = Threaded()))
    @test sim103.𝒪est[:l2] ≈ 5 atol = testTol
    sim105 = test_convergence(dts, prob, VCABM4(thread = Threaded()))
    @test sim105.𝒪est[:l2] ≈ 4 atol = testTol
end

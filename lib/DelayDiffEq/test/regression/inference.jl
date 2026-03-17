using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqRosenbrock
using Test

@testset "init" begin
    prob = prob_dde_constant_1delay_ip
    prob_scalar = prob_dde_constant_1delay_scalar

    inferred = [BS3(), Tsit5(), RK4(), Vern6()]
    for alg in inferred
        ddealg = MethodOfSteps(alg)

        @test_broken @inferred init(prob, ddealg)
        @test_broken @inferred init(prob_scalar, ddealg)
    end

    notinferred = [SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4()]
    for alg in notinferred
        ddealg = MethodOfSteps(alg)

        @test_broken @inferred init(prob, ddealg)
        @test_broken @inferred init(prob_scalar, ddealg)
    end
end

@testset "discontinuity_time" begin
    prob_inplace = prob_dde_constant_1delay_ip
    prob_scalar = prob_dde_constant_1delay_scalar

    for prob in (prob_inplace, prob_scalar)
        int = init(prob, MethodOfSteps(Tsit5()))
        @inferred DelayDiffEq.discontinuity_time(int, (u, p, t) -> 1.0, 0.0, (0.5, 1.5))
    end
end

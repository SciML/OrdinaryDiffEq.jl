using OrdinaryDiffEq, ODEProblemLibrary
using Test

@testset "init" begin
    prob = ODEProblemLibrary.prob_ode_linear
    prob2D = ODEProblemLibrary.prob_ode_2Dlinear

    inferred = [BS3(), Tsit5(), RK4(), Vern6()]
    for alg in inferred
        @inferred init(prob, alg)
        @inferred init(prob2D, alg)
    end

    #notinferred = [SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4()]
    #for alg in notinferred
    #    @test_broken @inferred init(prob, alg).t[1] == 0.0
    #    @test_broken @inferred init(prob2D, alg).t[1] == 0.0
    #end
end

using OrdinaryDiffEq, ODEProblemLibrary
using Test

@testset "init" begin
    prob = ODEProblemLibrary.prob_ode_linear
    prob2D = ODEProblemLibrary.prob_ode_2Dlinear

    inferred = [BS3(), Tsit5(), RK4(), Vern6(), SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4(), QNDF1(), QNDF()]
    for alg in inferred
        @inferred init(prob, alg)
        @inferred init(prob2D, alg)
    end
end

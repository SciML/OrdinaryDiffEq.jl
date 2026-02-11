using OrdinaryDiffEq, ODEProblemLibrary, ADTypes
using Test

@testset "init" begin
    prob = ODEProblemLibrary.prob_ode_linear
    prob2D = ODEProblemLibrary.prob_ode_2Dlinear

    inferred1 = [BS3(), Tsit5(), RK4(), Vern6()]
    for alg in inferred1
        @inferred init(prob, alg)
        @inferred init(prob2D, alg)
    end

    # ForwardDiff is not fully inferable
    autodiff = ADTypes.AutoFiniteDiff()
    inferred2 = [SDIRK2(; autodiff), TRBDF2(; autodiff), KenCarp4(; autodiff), Rosenbrock23(; autodiff), Rodas4(; autodiff)]
    for alg in inferred2
        @inferred init(prob, alg)
        @inferred init(prob2D, alg)
    end
end

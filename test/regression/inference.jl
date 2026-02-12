using OrdinaryDiffEq, ODEProblemLibrary, ADTypes
using Test

@testset "init" begin
    prob = ODEProblemLibrary.prob_ode_linear
    prob2D = ODEProblemLibrary.prob_ode_2Dlinear

    @testset "non-stiff" begin
        inferred1 = [BS3(), Tsit5(), RK4(), Vern6()]
        for alg in inferred1
            @inferred init(prob, alg)
            @inferred init(prob2D, alg)
        end
    end

    @testset "stiff default" begin
        # Stiff solvers are not fully inferable for the 2D problem with the default args
        inferred2 = [SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4()]
        for alg in inferred2
            @inferred init(prob, alg)
            @test_broken @inferred init(prob2D, alg)
        end
    end

    @testset "stiff fixed chunksize" begin
        # When choosing a fixed chunksize it works
        autodiff = ADTypes.AutoForwardDiff(; chunksize = 10)
        inferred3 = [SDIRK2(; autodiff), TRBDF2(; autodiff), KenCarp4(; autodiff), Rosenbrock23(; autodiff), Rodas4(; autodiff)]
        for alg in inferred3
            @inferred init(prob, alg)
            # In Julia < v1.12 only some of these are inferable
            if VERSION >= v"1.12"
                @inferred init(prob2D, alg)
            end
        end
    end

    @testset "stiff finite diff" begin
        # FiniteDiff works
        autodiff = ADTypes.AutoFiniteDiff()
        inferred4 = [SDIRK2(; autodiff), TRBDF2(; autodiff), KenCarp4(; autodiff), Rosenbrock23(; autodiff), Rodas4(; autodiff)]
        for alg in inferred4
            @inferred init(prob, alg)
            @inferred init(prob2D, alg)
        end
    end
end

using OrdinaryDiffEq, ODEProblemLibrary, ADTypes
using Test
using OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK

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
        inferred2 = [SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4()]
        for alg in inferred2
            @inferred init(prob, alg)
            @inferred init(prob2D, alg)
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

# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3200
# AutoVern7 + Rodas5P with both autodiff and linsolve should be inferrable
@testset "AutoVern7 + Rodas5P inference with autodiff and linsolve (#3200)" begin
    using StaticArrays, ADTypes, LinearSolve
    f(u, p, t) = SVector(-p.Ka * u[1], p.Ka * u[1] - p.CL * u[2] / p.Vc)
    u0 = SVector(0.0, 0.0)
    prob = ODEProblem(f, u0, (0.0, 10.0), (Ka = 1.0, CL = 1.0, Vc = 1.0))

    @inferred solve(
        prob,
        AutoVern7(
            Rodas5P(
                autodiff = AutoForwardDiff(chunksize = 1),
                linsolve = GenericLUFactorization()
            )
        )
    )
end

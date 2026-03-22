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
    #    @test_broken @inferred init(prob, alg)
    #    @test_broken @inferred init(prob2D, alg)
    #end
end

# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3200
# AutoVern7 + Rodas5P with both autodiff and linsolve should be inferable
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

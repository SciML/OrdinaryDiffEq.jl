using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

@testset "Constant delays" begin
    prob = DDEProblemLibrary.prob_dde_constant_2delays_ip
    prob_scalar = DDEProblemLibrary.prob_dde_constant_2delays_scalar

    algdict = Dict(BS3() => 2.4e-6,
                   Tsit5() => 4.5e-3,
                   RK4() => 1.1e-4,
                   Vern6() => 1.0e-3,
                   SDIRK2() => 2.3e-1,
                   TRBDF2() => 6.2e-2,
                   KenCarp4() => 5.6e-2,
                   Rosenbrock23() => 6.5e-4,
                   Rodas4() => 5.4e-4)

    for (alg, error) in algdict
        ddealg = MethodOfSteps(alg)

        sol = solve(prob, ddealg)
        @test sol.errors[:l∞] < error

        sol_scalar = solve(prob_scalar, ddealg)

        if alg isa OrdinaryDiffEq.OrdinaryDiffEqNewtonAdaptiveAlgorithm
            @test_broken sol.t ≈ sol_scalar.t
            @test_broken sol[1, :] ≈ sol_scalar.u
        else
            @test sol.t ≈ sol_scalar.t
            @test sol[1, :] ≈ sol_scalar.u
        end
    end
end
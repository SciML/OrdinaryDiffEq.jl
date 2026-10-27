using GlobalDiffEq, OrdinaryDiffEqTsit5
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqRosenbrock
using ADTypes: AutoForwardDiff
using Test
import DiffEqBase: SciMLBase

@testset "Algorithm traits forwarding" begin
    alg_inner = SSPRK33()
    alg = GlobalRichardson(alg_inner)

    @test SciMLBase.allows_arbitrary_number_types(alg) ==
        SciMLBase.allows_arbitrary_number_types(alg_inner)
    @test SciMLBase.allowscomplex(alg) == SciMLBase.allowscomplex(alg_inner)
    @test SciMLBase.isautodifferentiable(alg) == SciMLBase.isautodifferentiable(alg_inner)

    # The outer solve wraps f for these traits; an implicit inner solver needs the dual
    # signatures for its Jacobian.
    for inner in (SSPRK33(), Rodas5P(), Rodas5P(autodiff = AutoForwardDiff(chunksize = 2)))
        for alg in (GlobalRichardson(inner), GlobalErrorEstimation(inner))
            @test SciMLBase.forwarddiffs_model(alg) == SciMLBase.forwarddiffs_model(inner)
            @test SciMLBase.forwarddiff_chunksize(alg) == SciMLBase.forwarddiff_chunksize(inner)
        end
    end
    @test SciMLBase.forwarddiffs_model(GlobalErrorEstimation(Rodas5P()))
    decay!(du, u, p, t) = (du[1] = -u[1]; nothing)
    prob = ODEProblem(decay!, [1.0], (0.0, 1.0))
    @test SciMLBase.successful_retcode(solve(prob, GlobalRichardson(Rodas5P()); dt = 0.01))
end

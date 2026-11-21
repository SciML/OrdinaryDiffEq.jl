using GlobalDiffEq, OrdinaryDiffEqTsit5
using OrdinaryDiffEqSSPRK
using Test
import DiffEqBase: SciMLBase

@testset "Algorithm traits forwarding" begin
    alg_inner = SSPRK33()
    alg = GlobalRichardson(alg_inner)

    @test SciMLBase.allows_arbitrary_number_types(alg) ==
        SciMLBase.allows_arbitrary_number_types(alg_inner)
    @test SciMLBase.allowscomplex(alg) == SciMLBase.allowscomplex(alg_inner)
    @test SciMLBase.isautodifferentiable(alg) == SciMLBase.isautodifferentiable(alg_inner)
end

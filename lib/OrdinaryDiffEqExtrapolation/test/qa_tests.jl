using OrdinaryDiffEqExtrapolation
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqExtrapolation
    )
end
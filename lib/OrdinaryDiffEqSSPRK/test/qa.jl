using OrdinaryDiffEqSSPRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqSSPRK
    )
end
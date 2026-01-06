using OrdinaryDiffEqLowOrderRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqLowOrderRK
    )
end

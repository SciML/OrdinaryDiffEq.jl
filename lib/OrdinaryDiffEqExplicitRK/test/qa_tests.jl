using OrdinaryDiffEqExplicitRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqExplicitRK
    )
end
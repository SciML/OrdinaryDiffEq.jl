using OrdinaryDiffEqStabilizedRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqStabilizedRK
    )
end

using OrdinaryDiffEqStabilizedIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqStabilizedIRK
    )
end

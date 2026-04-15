using OrdinaryDiffEqPDIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqPDIRK
    )
end

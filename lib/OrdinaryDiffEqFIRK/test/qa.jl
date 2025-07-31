using OrdinaryDiffEqFIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFIRK
    )
end

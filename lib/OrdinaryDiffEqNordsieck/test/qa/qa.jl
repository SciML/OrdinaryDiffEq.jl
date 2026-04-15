using OrdinaryDiffEqNordsieck
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqNordsieck
    )
end

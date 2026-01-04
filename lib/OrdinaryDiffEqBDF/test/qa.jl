using OrdinaryDiffEqBDF
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqBDF
    )
end

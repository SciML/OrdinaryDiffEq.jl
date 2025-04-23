using OrdinaryDiffEqCore
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqCore
    )
end
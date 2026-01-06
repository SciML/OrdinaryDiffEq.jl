using OrdinaryDiffEqLinear
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqLinear
    )
end

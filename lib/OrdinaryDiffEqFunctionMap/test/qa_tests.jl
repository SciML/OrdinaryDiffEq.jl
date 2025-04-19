using OrdinaryDiffEqFunctionMap
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFunctionMap
    )
end
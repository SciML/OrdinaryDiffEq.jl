using OrdinaryDiffEqFunctionMap
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFunctionMap;
        piracies = false  # Piracy is necessary for default algorithm dispatch
    )
end

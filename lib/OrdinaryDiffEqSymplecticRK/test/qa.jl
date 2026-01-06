using OrdinaryDiffEqSymplecticRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqSymplecticRK
    )
end

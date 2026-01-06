using OrdinaryDiffEqRKN
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqRKN
    )
end

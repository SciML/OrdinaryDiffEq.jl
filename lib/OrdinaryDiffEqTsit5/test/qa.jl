using OrdinaryDiffEqTsit5
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqTsit5
    )
end

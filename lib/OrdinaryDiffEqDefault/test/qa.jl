using OrdinaryDiffEqDefault
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqDefault;
        piracies = false
    )
end
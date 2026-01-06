using OrdinaryDiffEqFeagin
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFeagin
    )
end

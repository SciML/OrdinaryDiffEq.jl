using OrdinaryDiffEqQPRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqQPRK
    )
end

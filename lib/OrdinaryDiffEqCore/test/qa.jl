using OrdinaryDiffEqCore
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqCore;
        piracies = false,
        unbound_args = false
    )
end
using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqCore
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqCore;
        piracies = false,
        unbound_args = false
    )
end

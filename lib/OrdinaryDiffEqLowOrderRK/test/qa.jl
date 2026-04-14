using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqLowOrderRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqLowOrderRK
    )
end

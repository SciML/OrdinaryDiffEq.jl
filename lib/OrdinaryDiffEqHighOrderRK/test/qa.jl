using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqHighOrderRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqHighOrderRK
    )
end

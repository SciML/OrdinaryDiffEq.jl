using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqVerner
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqVerner
    )
end

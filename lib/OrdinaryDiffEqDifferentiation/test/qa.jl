using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqDifferentiation
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqDifferentiation;
        piracies = false,
        ambiguities = false
    )
end

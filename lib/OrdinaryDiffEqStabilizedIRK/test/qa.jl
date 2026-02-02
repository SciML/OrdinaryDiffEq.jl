using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqStabilizedIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqStabilizedIRK
    )
end

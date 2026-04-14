using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqFIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFIRK
    )
end

using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqSDIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqSDIRK
    )
end

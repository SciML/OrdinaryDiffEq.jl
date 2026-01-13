using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqExponentialRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqExponentialRK
    )
end

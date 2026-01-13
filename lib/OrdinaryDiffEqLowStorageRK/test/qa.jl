using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqLowStorageRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqLowStorageRK
    )
end

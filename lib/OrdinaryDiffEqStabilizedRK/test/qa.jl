using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqStabilizedRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqStabilizedRK;
        deps_compat = (check_extras = false,)
    )
end

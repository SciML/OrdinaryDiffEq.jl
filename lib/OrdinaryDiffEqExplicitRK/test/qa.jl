using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqExplicitRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqExplicitRK;
        deps_compat = (check_extras = false,)
    )
end

using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqHighOrderRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqHighOrderRK;
        deps_compat = (check_extras = false,)
    )
end

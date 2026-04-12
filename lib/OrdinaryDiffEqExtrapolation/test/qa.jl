using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqExtrapolation
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqExtrapolation;
        deps_compat = (check_extras = false,)
    )
end

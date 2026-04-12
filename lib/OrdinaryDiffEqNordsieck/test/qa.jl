using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqNordsieck
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqNordsieck;
        deps_compat = (check_extras = false,)
    )
end

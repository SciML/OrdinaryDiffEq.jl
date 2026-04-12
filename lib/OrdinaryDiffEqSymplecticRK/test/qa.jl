using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqSymplecticRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqSymplecticRK;
        deps_compat = (check_extras = false,)
    )
end

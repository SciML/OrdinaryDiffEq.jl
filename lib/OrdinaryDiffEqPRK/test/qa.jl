using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqPRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqPRK;
        deps_compat = (check_extras = false,)
    )
end

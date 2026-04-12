using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqRKN
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqRKN;
        deps_compat = (check_extras = false,)
    )
end

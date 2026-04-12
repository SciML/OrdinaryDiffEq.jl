using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqFeagin
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqFeagin;
        deps_compat = (check_extras = false,)
    )
end

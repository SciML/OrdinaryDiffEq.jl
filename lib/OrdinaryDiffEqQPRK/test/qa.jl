using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqQPRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqQPRK;
        deps_compat = (check_extras = false,)
    )
end

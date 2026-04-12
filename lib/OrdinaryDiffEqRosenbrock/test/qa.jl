using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqRosenbrock
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqRosenbrock;
        deps_compat = (check_extras = false,)
    )
end

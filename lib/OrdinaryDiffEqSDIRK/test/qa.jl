using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqSDIRK
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqSDIRK;
        deps_compat = (check_extras = false,)
    )
end

using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqAdamsBashforthMoulton
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqAdamsBashforthMoulton;
        deps_compat = (check_extras = false,)
    )
end

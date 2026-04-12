using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqIMEXMultistep
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqIMEXMultistep;
        deps_compat = (check_extras = false,)
    )
end

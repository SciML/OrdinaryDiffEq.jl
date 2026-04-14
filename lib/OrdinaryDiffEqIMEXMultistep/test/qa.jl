using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqIMEXMultistep
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqIMEXMultistep
    )
end

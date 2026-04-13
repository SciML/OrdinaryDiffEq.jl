using Pkg
Pkg.add("Aqua")

using ImplicitDiscreteSolve
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        ImplicitDiscreteSolve;
        piracies = false
    )
end

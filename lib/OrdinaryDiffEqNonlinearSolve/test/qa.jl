using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqNonlinearSolve
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqNonlinearSolve;
        piracies = false,
        deps_compat = (check_extras = false,)
    )
end

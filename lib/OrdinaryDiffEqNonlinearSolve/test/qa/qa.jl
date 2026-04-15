using OrdinaryDiffEqNonlinearSolve
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqNonlinearSolve;
        piracies = false
    )
end

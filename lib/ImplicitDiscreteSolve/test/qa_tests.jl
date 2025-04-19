using ImplicitDiscreteSolve
using Aqua


@testset "Aqua" begin
    Aqua.test_all(
        ImplicitDiscreteSolve
    )
end
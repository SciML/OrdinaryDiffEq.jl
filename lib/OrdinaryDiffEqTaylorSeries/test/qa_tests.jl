using OrdinaryDiffEqTaylorSeries
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqTaylorSeries
    )
end
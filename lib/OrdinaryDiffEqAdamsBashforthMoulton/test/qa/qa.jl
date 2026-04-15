using OrdinaryDiffEqAdamsBashforthMoulton
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqAdamsBashforthMoulton
    )
end

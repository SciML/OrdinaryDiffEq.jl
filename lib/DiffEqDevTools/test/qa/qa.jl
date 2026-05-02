using DiffEqDevTools
using Aqua
using Test

@testset "Aqua" begin
    Aqua.test_all(
        DiffEqDevTools;
        ambiguities = false,
        piracies = false,
        unbound_args = false,
        stale_deps = false,
        deps_compat = (; check_extras = false),
    )
end

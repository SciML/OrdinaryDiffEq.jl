using OrdinaryDiffEqTaylorSeries
using Aqua

@testset "Aqua" begin
    Aqua.test_all(
        OrdinaryDiffEqTaylorSeries;
        unbound_args = false,
        undefined_exports = false,
        stale_deps = false,
        deps_compat = false,
        ambiguities = false
    )
end

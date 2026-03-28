using DelayDiffEq, Test
import SciMLBase

@testset "Aqua Tests" begin
    using Aqua

    Aqua.test_all(
        DelayDiffEq; ambiguities = false, piracies = false,
        stale_deps = false, deps_compat = false
    )
    Aqua.test_ambiguities(DelayDiffEq; recursive = false)
    Aqua.test_stale_deps(DelayDiffEq)
    Aqua.test_deps_compat(DelayDiffEq)
    # Allow piracy for the default solver methods and SDE integration
    using OrdinaryDiffEqCore: StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        StochasticDiffEqConstantCache, StochasticDiffEqMutableCache
    Aqua.test_piracies(
        DelayDiffEq;
        treat_as_own = [
            SciMLBase.DDEProblem,
            StochasticDiffEqAlgorithm,
            StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqConstantCache,
            StochasticDiffEqMutableCache,
        ]
    )
end

@testset "Explicit Imports Tests" begin
    using ExplicitImports

    @test check_no_implicit_imports(DelayDiffEq; skip = (Base, Core), ignore = (Symbol("@reexport"),)) ===
        nothing
    @test check_no_stale_explicit_imports(DelayDiffEq, ignore = (:AbstractVerbositySpecifier, :Standard)) === nothing
    @test check_all_qualified_accesses_via_owners(DelayDiffEq) === nothing
end

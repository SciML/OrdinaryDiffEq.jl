using ExplicitImports, OrdinaryDiffEq
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(
        OrdinaryDiffEq; skip = (Base, Core, SciMLBase)
    ) === nothing

    @test check_no_stale_explicit_imports(OrdinaryDiffEq) === nothing

    @test check_all_qualified_accesses_via_owners(OrdinaryDiffEq) === nothing
end

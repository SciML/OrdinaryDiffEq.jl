using SciMLTesting, GlobalDiffEq, Test
using JET
using OrdinaryDiffEqTsit5, OrdinaryDiffEqSSPRK

run_qa(
    GlobalDiffEq;
    reexports_allow = union(public_api_names(DiffEqBase), (:DiffEqBase,)),
    explicit_imports = true,
    # GlobalDiffEq's rendered documentation lives in the monorepo docs, two
    # directories up from the sublibrary root.
    api_docs_kwargs = (;
        docs_src = joinpath(dirname(dirname(pkgdir(GlobalDiffEq))), "docs", "src"),
    ),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            # `SciMLBase.__solve` is SciMLBase's internal solve entry point (not
            # part of the public API); GlobalDiffEq overloads it via its owner
            # SciMLBase.
            ignore = (:__solve,),
        ),
    ),
    # `@reexport using DiffEqBase` deliberately reexports DiffEqBase's API, so
    # `DiffEqBase`, `ODEProblem`, and `solve` are inherently implicit. Tracked in
    # https://github.com/SciML/GlobalDiffEq.jl/issues/53
    ei_broken = (:no_implicit_imports,),
)

@testset "GlobalRichardson static analysis" begin
    @test_opt GlobalRichardson(SSPRK33())
    @test GlobalRichardson{typeof(SSPRK33())} <: GlobalDiffEq.GlobalDiffEqAlgorithm
    @test GlobalRichardson(Tsit5()) isa GlobalDiffEq.GlobalDiffEqAlgorithm
end

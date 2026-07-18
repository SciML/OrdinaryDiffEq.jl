using SciMLTesting, GlobalDiffEq, Test
using JET
using OrdinaryDiffEqTsit5, OrdinaryDiffEqSSPRK

# `@reexport using DiffEqBase` republishes DiffEqBase's API; those names are
# documented and rendered at their owning packages, not in the OrdinaryDiffEq
# manual's GlobalDiffEq page, so only GlobalDiffEq's own names are held to the
# rendered-docs check.
reexported_names = Tuple(
    filter(public_api_names(GlobalDiffEq)) do name
        object = getfield(GlobalDiffEq, name)
        !(object isa Union{Function, Type, Module}) ||
            parentmodule(object) !== GlobalDiffEq
    end
)

run_qa(
    GlobalDiffEq;
    explicit_imports = true,
    # GlobalDiffEq's rendered documentation lives in the monorepo docs, two
    # directories up from the sublibrary root.
    api_docs_kwargs = (;
        rendered = true,
        docs_src = joinpath(dirname(dirname(pkgdir(GlobalDiffEq))), "docs", "src"),
        rendered_ignore = reexported_names,
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

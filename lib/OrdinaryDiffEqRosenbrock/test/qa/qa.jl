using SciMLTesting, OrdinaryDiffEqRosenbrock, Test

# OrdinaryDiffEqRosenbrock consumes the solver-author extension API of its sibling
# packages, most of which is now declared public (so those names are no longer
# listed here). The names below are the genuine residual: cross-package internals
# whose *owner* does not make them public, plus a couple of owner-internal Core
# helpers. Each is grouped and commented by owning package so that any name later
# made public will surface as an Unexpected Pass and can be dropped. See
# SciML/OrdinaryDiffEq.jl#3776.

# Explicitly-imported names that are non-public in (and/or owned outside) the
# module they are imported from.
const ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS = (
    :_unwrap_val, :_vec,                     # owner SciMLBase, non-public internals
    :calculate_residuals, :calculate_residuals!, :initialize!,  # owner DiffEqBase, non-public
    :copyat_or_push!,                        # owner RecursiveArrayTools, non-public
    :trivial_limiter!,                       # OrdinaryDiffEqCore owner-internal, non-public
)

# Qualified accesses to names that are non-public in (and/or owned outside) the
# module they are accessed through.
const ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES = (
    Symbol("@set"),        # owner Accessors, accessed via SciMLBase.@set, non-public
    :Cartesian,            # Base's bounded unroll macro is required by the fused kernels.
    :setindex,             # Base internal, non-public
)

run_qa(
    OrdinaryDiffEqRosenbrock;
    api_docs_kwargs = (; docs_src = joinpath(pkgdir(OrdinaryDiffEqRosenbrock), "..", "..", "docs", "src")),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (
        # `@reexport using SciMLBase` makes this explicitly imported public macro
        # appear implicit to ExplicitImports.
        no_implicit_imports = (; ignore = (Symbol("@def"),)),
        all_explicit_imports_via_owners = (; ignore = ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS),
        all_explicit_imports_are_public = (; ignore = ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS),
        all_qualified_accesses_via_owners = (; ignore = ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES),
        all_qualified_accesses_are_public = (; ignore = ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES),
    ),
)

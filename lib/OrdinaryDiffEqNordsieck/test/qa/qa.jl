using SciMLTesting, OrdinaryDiffEqNordsieck, Test

run_qa(
    OrdinaryDiffEqNordsieck;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` is the package's public-API surface; the
        # reexport idiom is an unavoidable implicit import of these three names.
        no_implicit_imports = (;
            ignore = (:SciMLBase, :Reexport, Symbol("@reexport")),
        ),
        # Genuine solver-author internals that their owners deliberately keep
        # non-public. Everything OrdinaryDiffEqCore made public for the solver
        # extension API has been dropped; drop each remaining entry once its
        # owner marks it public upstream (tracked in SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals (not in the public extension API)
                :_resolved_QT, :resolve_basic,
                :get_current_alg_order,
                :qsteady_min_default, :qsteady_max_default,
                :ode_interpolant, :ode_interpolant!,
                :trivial_limiter!,
                # OrdinaryDiffEqTsit5 first-step caches (owner-internal)
                :Tsit5Cache, :Tsit5ConstantCache,
            ),
        ),
    ),
)

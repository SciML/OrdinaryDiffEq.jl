using SciMLTesting, OrdinaryDiffEqVerner, Test

run_qa(
    OrdinaryDiffEqVerner;
    explicit_imports = true,
    ei_kwargs = (
        # `SciMLBase` module name is brought into scope by `@reexport using SciMLBase`
        # (intentional: this package re-exports the SciMLBase public API and also uses
        # it as a qualified accessor, e.g. `SciMLBase.solve`, `SciMLBase.AutoSpecialize`).
        no_implicit_imports = (; ignore = (:SciMLBase,)),
        # Qualified accesses to names that are not (yet) declared public by their owners.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase internals (not public on registered SciMLBase 3.30.x)
                :has_lazy_interpolation,
                # OrdinaryDiffEqCore test-only helpers (deliberately kept non-public)
                :lorenz, :lorenz_oop,
            )),
        # Explicit imports of names not (yet) declared public by their owners.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal names (private codegen macros +
                # cache/interpolation internals deliberately not part of the public
                # solver-author API declared in PHASE A)
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"),
                :CompiledFloats, :DerivativeOrderNotPossibleError,
                :_ode_interpolant!, :trivial_limiter!,
                # SciMLBase internals (not public on registered SciMLBase 3.30.x)
                Symbol("@def"), :_unwrap_val,
                # TruncatedStacktraces internal macro
                Symbol("@truncate_stacktrace"),
            )),
    ),
)

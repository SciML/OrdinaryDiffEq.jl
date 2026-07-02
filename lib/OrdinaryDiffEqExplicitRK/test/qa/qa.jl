using SciMLTesting, OrdinaryDiffEqExplicitRK, Test

run_qa(
    OrdinaryDiffEqExplicitRK;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore solver-internal names still kept non-public
                # (the `!`/error variants were not part of the make-public surface).
                :DerivativeOrderNotPossibleError, :_ode_interpolant!,
                # TruncatedStacktraces-owned macro (external, not `public`).
                Symbol("@truncate_stacktrace"),
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore solver-internal interpolation hooks kept non-public.
                :interpolation_differential_vars, :hermite_interpolant!,
                # Base.Cartesian codegen macros (external, not `public`).
                Symbol("@nexprs"), Symbol("@nif"),
            ),
        ),
    ),
)

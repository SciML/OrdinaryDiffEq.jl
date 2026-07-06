using SciMLTesting, OrdinaryDiffEqLowOrderRK, Test

run_qa(
    OrdinaryDiffEqLowOrderRK;
    explicit_imports = true,
    ei_kwargs = (;
        # SciMLBase's has_lazy_interpolation interpolation trait, accessed via its
        # owner module. Not declared public on the registered SciMLBase release.
        all_qualified_accesses_are_public = (;
            ignore = (
                :has_lazy_interpolation,  # SciMLBase
            ),
        ),
        # Owner-internal names imported from their owning module but deliberately
        # NOT declared `public`. OrdinaryDiffEqCore's public solver-author API is
        # now declared public on this branch; the residual here is the codegen/perf
        # macros and a few helpers/types that OrdinaryDiffEqCore keeps internal,
        # DiffEqBase's private @tight_loop_macros, and SciMLBase internals.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore (owner-internal; not part of the public
                # solver-author surface declared in OrdinaryDiffEqCore)
                :CompiledFloats, :DerivativeOrderNotPossibleError,
                Symbol("@fold"), :_ode_interpolant!,
                Symbol("@OnDemandTableauExtract"),
                :ssp_coefficient, :trivial_limiter!,
                # DiffEqBase (private loop macro)
                Symbol("@tight_loop_macros"),
                # SciMLBase
                Symbol("@def"), :_unwrap_val,
            ),
        ),
    ),
)

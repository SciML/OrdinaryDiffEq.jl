using SciMLTesting, StochasticDiffEqImplicit, Test
using JET

run_qa(
    StochasticDiffEqImplicit;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` is owned by FastBroadcast and re-exported through DiffEqBase; FastBroadcast
        # is not a direct dependency, so it is imported from DiffEqBase by design.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqNonlinearSolve — owner-internal, no public alternative
                :build_nlsolver, :markfirststage!, :nlsolve!, :nlsolvefail,
                # `@..` (FastBroadcast macro re-exported via DiffEqBase) is not `public` there.
                Symbol("@.."),
                # StochasticDiffEqCore internal codegen macro (owner-internal).
                Symbol("@cache"),
                # OrdinaryDiffEqCore internals with no public replacement yet.
                :current_extrapolant!, :isnewton, :_fixup_ad,
                # SciMLBase internals (not yet `public` on registered SciMLBase).
                :has_Wfact, :_vec, :_reshape, :_unwrap_val,
            ),
        ),
    ),
)

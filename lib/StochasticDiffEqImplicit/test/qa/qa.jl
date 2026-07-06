using SciMLTesting, StochasticDiffEqImplicit, Test
using JET

run_qa(
    StochasticDiffEqImplicit;
    # Scope JET to this package in `:typo` mode, matching the OrdinaryDiffEq* solver
    # sublibraries. The deprecated `target_defined_modules = true` also reported
    # `nlsolve!`/`get_W`/`set_new_W!` "no matching method" cross-package false
    # positives that arise only because `perform_step!` leaves `integrator`
    # untyped (at runtime it is always an `ODEIntegrator <: DEIntegrator`, so the
    # functional Core tests exercise these paths without error).
    jet_kwargs = (; target_modules = (StochasticDiffEqImplicit,), mode = :typo),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` is owned by FastBroadcast and re-exported through DiffEqBase; FastBroadcast
        # is not a direct dependency, so it is imported from DiffEqBase by design.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        all_explicit_imports_are_public = (;
            ignore = (
                # `@..` (FastBroadcast macro re-exported via DiffEqBase) is not `public` there.
                Symbol("@.."),
                # StochasticDiffEqCore internal codegen macro (owner-internal).
                Symbol("@cache"),
                # OrdinaryDiffEqCore internals with no public replacement yet.
                :current_extrapolant!, :isnewton, :_fixup_ad,
                # SciMLBase internals (not yet `public` on registered SciMLBase).
                :has_Wfact, :_vec, :_reshape, :_unwrap_val,
            )),
    ),
)

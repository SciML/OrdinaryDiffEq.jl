using SciMLTesting, OrdinaryDiffEqTaylorSeries, Test

# Every name in the ignore lists below is genuinely internal (not exported and
# not declared `public`) in its OWNER package, so there is no public name to
# migrate to. The base solver-author API (perform_step!/alg_cache/controllers/
# etc.) is now declared `public` in OrdinaryDiffEqCore, so those entries have
# been dropped; only the genuine residual remains (see SciML/OrdinaryDiffEq.jl#3776).
run_qa(
    OrdinaryDiffEqTaylorSeries;
    aqua_kwargs = (;
        unbound_args = false, undefined_exports = false, stale_deps = false,
        ambiguities = false,
    ),
    explicit_imports = true,
    ei_kwargs = (;
        # Explicit imports of names that are not `public` in their owner.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore limiter default, deliberately kept
                # non-public (codegen/perf helper, not extension API).
                :trivial_limiter!,
                # other-package internals with no public alias
                Symbol("@truncate_stacktrace"), :FunctionWrapper,
            ),
        ),
        # Qualified accesses (Module.name) of non-public names.
        all_qualified_accesses_are_public = (;
            ignore = (
                # lorenz/lorenz_oop precompile-workload helpers (OrdinaryDiffEqCore,
                # deliberately non-public per its public-API comment).
                :lorenz, :lorenz_oop,
                # TaylorDiff internals (no public API for these primitives)
                :append_coefficient, :flatten, :get_coefficient, :make_seed,
                :partials,
                # Base / Symbolics internals
                :RefValue, :variables,
            ),
        ),
    ),
)

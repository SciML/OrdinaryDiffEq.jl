using SciMLTesting, OrdinaryDiffEqTaylorSeries, Test

# Every name in the ignore lists below is genuinely internal (not exported and
# not declared `public`) in its OWNER package, so there is no public name to
# migrate to. The vast majority are OrdinaryDiffEqCore solver-interface hooks
# (perform_step!/alg_cache/initialize!/controllers/etc.) that sublibraries must
# extend; once OrdinaryDiffEqCore marks them `public` these entries should be
# dropped (see SciML/OrdinaryDiffEq.jl#3776).
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
                # OrdinaryDiffEqCore solver-interface internals
                Symbol("@cache"), :DerivativeOrderNotPossibleError,
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!,
                :alg_cache, :alg_stability_size, :explicit_rk_docstring,
                :get_current_adaptive_order, :get_current_alg_order,
                :get_fsalfirstlast, :perform_step!, :step_accept_controller!,
                :stepsize_controller!, :trivial_limiter!, :unwrap_alg,
                # DiffEqBase residual-calculation internals
                :calculate_residuals, :calculate_residuals!,
                # other-package internals with no public alias
                Symbol("@truncate_stacktrace"), :FunctionWrapper,
            ),
        ),
        # Qualified accesses (Module.name) of non-public names.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore stats/controller internals
                :get_EEst, :increment_nf!, :set_EEst!, :sync_controllers!,
                # lorenz/lorenz_oop test functions used in @compile_workload
                :lorenz, :lorenz_oop,
                # TaylorDiff internals (no public API for these primitives)
                :append_coefficient, :flatten, :get_coefficient, :make_seed,
                :partials,
                # SciMLBase / Base / Symbolics internals
                :NoSpecialize, :RefValue, :variables,
            ),
        ),
    ),
)

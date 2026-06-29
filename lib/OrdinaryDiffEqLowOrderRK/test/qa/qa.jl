using SciMLTesting, OrdinaryDiffEqLowOrderRK, Test

run_qa(
    OrdinaryDiffEqLowOrderRK;
    explicit_imports = true,
    ei_kwargs = (;
        # Solver-internal names with no public alternative: OrdinaryDiffEqCore's
        # increment_nf!/set_EEst! step-statistics mutators and SciMLBase's
        # has_lazy_interpolation interpolation trait. Accessed via qualified
        # owner module; none are public on the registered releases.
        all_qualified_accesses_are_public = (;
            ignore = (
                :has_lazy_interpolation,  # SciMLBase
                :increment_nf!,           # OrdinaryDiffEqCore
                :set_EEst!,               # OrdinaryDiffEqCore
            ),
        ),
        # Core solver-interface names that are imported from their owner module
        # but not yet declared `public`. These are the OrdinaryDiffEq internal
        # extension/dispatch API (OrdinaryDiffEqCore), plus DiffEqBase residual
        # helpers and SciMLBase internals (@def/_unwrap_val). No public
        # alternative exists; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore
                :accept_step_controller, :alg_cache, :alg_stability_size,
                :AutoAlgSwitch, :beta1_default, :beta2_default,
                Symbol("@cache"), :CompiledFloats, :CompositeAlgorithm,
                :constvalue, :DerivativeOrderNotPossibleError,
                :explicit_rk_docstring, Symbol("@fold"),
                :generic_solver_docstring, :get_fsalfirstlast, :isfsal,
                :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!,
                Symbol("@OnDemandTableauExtract"),
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqAlgorithm,
                :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqExponentialAlgorithm,
                :OrdinaryDiffEqMutableCache, :perform_step!, :ssp_coefficient,
                :trivial_limiter!, :unwrap_alg,
                # DiffEqBase
                :calculate_residuals, :calculate_residuals!, :prepare_alg,
                Symbol("@tight_loop_macros"),
                # SciMLBase
                Symbol("@def"), :_unwrap_val,
            ),
        ),
    ),
)

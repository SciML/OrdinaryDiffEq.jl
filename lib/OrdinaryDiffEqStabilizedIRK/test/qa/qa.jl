using SciMLTesting, OrdinaryDiffEqStabilizedIRK, Test

run_qa(
    OrdinaryDiffEqStabilizedIRK;
    explicit_imports = true,
    ei_kwargs = (;
        # Solver-interface names imported from their owner module but not yet
        # declared `public` on the registered releases. These are the
        # OrdinaryDiffEq internal extension/dispatch API (OrdinaryDiffEqCore,
        # OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqDifferentiation), plus
        # DiffEqBase residual helpers and SciMLBase internals. No public
        # alternative exists; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore
                :default_controller, :IController, :gamma_default, :issplit,
                :perform_step!, :unwrap_alg, :fac_default_gamma,
                :OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                :OrdinaryDiffEqMutableCache, :OrdinaryDiffEqConstantCache,
                :alg_cache, Symbol("@cache"), :get_fsalfirstlast,
                :increment_nf!, :set_EEst!, :generic_solver_docstring,
                :_fixup_ad, :get_W, :isnewton,
                # OrdinaryDiffEqNonlinearSolve
                :NLNewton, :nlsolve!, :build_nlsolver, :markfirststage!,
                :du_alias_or_new,
                # OrdinaryDiffEqDifferentiation
                :dolinsolve, :update_W!,
                # DiffEqBase
                :calculate_residuals, :calculate_residuals!,
                # SciMLBase
                :_unwrap_val, :_reshape, :_vec,
            ),
        ),
    ),
)

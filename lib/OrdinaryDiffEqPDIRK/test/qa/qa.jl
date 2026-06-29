using SciMLTesting, OrdinaryDiffEqPDIRK, Test

run_qa(
    OrdinaryDiffEqPDIRK;
    explicit_imports = true,
    ei_kwargs = (;
        # Solver-internal names extended/dispatched from their owner modules but
        # not (yet) declared `public` on the registered releases. No public
        # alternative exists; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore
                Symbol("@cache"), Symbol("@threaded"),
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqNewtonAlgorithm, :_fixup_ad, :alg_cache,
                :constvalue, :differentiation_rk_docstring, :get_fsalfirstlast,
                :isfsal, :isthreaded, :perform_step!, :unwrap_alg,
                # OrdinaryDiffEqNonlinearSolve
                :NLNewton, :build_nlsolver, :markfirststage!, :nlsolve!,
                :nlsolvefail,
                # SciMLBase
                :_unwrap_val,
            ),
        ),
    ),
)

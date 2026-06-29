using SciMLTesting, OrdinaryDiffEqSDIRK, Test

# The remaining ignores are all names that are genuinely non-public internals of
# the SDIRK solver stack (OrdinaryDiffEqCore / OrdinaryDiffEqNonlinearSolve /
# OrdinaryDiffEqDifferentiation), plus a handful of non-public SciMLBase /
# DiffEqBase / ConstructionBase names with no public alternative. They are owned
# by the module they are imported/accessed from (so the *_via_owners checks pass);
# they are only flagged because those owners have not yet declared them `public`.
# Each will drop off once the owning package marks the corresponding name public.
run_qa(
    OrdinaryDiffEqSDIRK;
    explicit_imports = true,
    ei_kwargs = (
        # `constructorof` is owned by ConstructionBase (not a direct dep) and is
        # reached through SciMLBase's re-export; it is the only non-owner access.
        all_qualified_accesses_via_owners = (; ignore = (:constructorof,)),
        all_qualified_accesses_are_public = (;
            ignore = (
                # non-public in SciMLBase
                :constructorof, :NoSpecialize, :FunctionWrapperSpecialize,
                # non-public OrdinaryDiffEqCore internals
                :increment_nf!, :set_EEst!, :get_EEst, :strip_cache,
                :lorenz, :lorenz_oop,
            )),
        all_explicit_imports_are_public = (;
            ignore = (
                # non-public SciMLBase internals
                :_reshape, :_unwrap_val, :_vec,
                # non-public DiffEqBase internals
                :calculate_residuals, :calculate_residuals!,
                # non-public TruncatedStacktraces macro
                Symbol("@truncate_stacktrace"),
                # non-public OrdinaryDiffEqCore internals
                :CompiledFloats, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                :OrdinaryDiffEqNewtonAlgorithm,
                :_fixup_ad, :alg_cache, :alg_extrapolates, :constvalue,
                :current_extrapolant!, :generic_solver_docstring,
                :get_fsalfirstlast, :isesdirk, :issplit, :perform_step!,
                :ssp_coefficient, :trivial_limiter!, :unwrap_alg,
                :COEFFICIENT_MULTISTEP, :get_W, :isnewton, :set_new_W!,
                # non-public OrdinaryDiffEqNonlinearSolve internals
                :NLNewton, :build_nlsolver, :du_alias_or_new, :markfirststage!,
                :nlsolve!, :nlsolvefail,
                # non-public OrdinaryDiffEqDifferentiation internals
                :dolinsolve,
            )),
    ),
)

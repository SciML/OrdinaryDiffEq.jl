using SciMLTesting, OrdinaryDiffEqRKN, Test

# The only remaining EI exceptions are genuinely internal (non-public) names from
# the SciML solver stack that OrdinaryDiffEqRKN must extend/use but which have no
# public alias yet. They are owned by the modules they are imported/accessed from;
# the sole remaining issue is that those owners have not marked them `public`.
# See SciML/OrdinaryDiffEq.jl#3776 for the make-public tracking.
run_qa(
    OrdinaryDiffEqRKN;
    explicit_imports = true,
    ei_kwargs = (;
        # Qualified accesses of OrdinaryDiffEqCore internals (no public alias).
        all_qualified_accesses_are_public = (;
            ignore = (:increment_nf!, :set_EEst!),),
        # Explicit imports of names that are owned by the importing module but
        # not yet declared `public` there: OrdinaryDiffEqCore solver-extension
        # internals, the DiffEqBase residual interface + @tight_loop_macros, and
        # SciMLBase's @def.
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@cache"), :CompiledFloats, :OrdinaryDiffEqAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqAdaptivePartitionedAlgorithm,
                :OrdinaryDiffEqPartitionedAlgorithm, :_ode_interpolant,
                :_ode_interpolant!, :alg_cache, :alg_extrapolates, :constvalue,
                :generic_solver_docstring, :get_fsalfirstlast, :perform_step!,
                Symbol("@tight_loop_macros"), :calculate_residuals,
                :calculate_residuals!, Symbol("@def"),
            ),),
    ),
)

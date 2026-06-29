using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    explicit_imports = true,
    ei_kwargs = (;
        # OrdinaryDiffEqCore-internal (non-public) names the solver genuinely
        # needs from Core; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                :perform_step!,
                :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqPartitionedAlgorithm,
                :CompiledFloats,
                :alg_cache,
                Symbol("@cache"),
                :constvalue,
                :get_fsalfirstlast,
                :generic_solver_docstring,
                :default_linear_interpolation,
            ),
        ),
        # `increment_nf!` is accessed as `OrdinaryDiffEqCore.increment_nf!`;
        # non-public in Core, tracked for make-public in #3776.
        all_qualified_accesses_are_public = (; ignore = (:increment_nf!,)),
    ),
)

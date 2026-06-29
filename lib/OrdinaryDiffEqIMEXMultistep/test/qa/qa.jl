using SciMLTesting, OrdinaryDiffEqIMEXMultistep, Test

run_qa(
    OrdinaryDiffEqIMEXMultistep;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` brings the `SciMLBase` module name into scope
        # implicitly; this is the intended SciML re-export pattern, not a stray import.
        no_implicit_imports = (; ignore = (:SciMLBase,)),
        # `OrdinaryDiffEqCore.increment_nf!` is a solver-internal stats helper that is
        # not (yet) declared public in OrdinaryDiffEqCore.
        all_qualified_accesses_are_public = (; ignore = (:increment_nf!,)),
        # Solver-internal building blocks imported from their owning packages but not
        # (yet) declared public there: OrdinaryDiffEqCore cache/algorithm/dispatch
        # internals, OrdinaryDiffEqNonlinearSolve nlsolver internals, and the
        # SciMLBase-internal `_unwrap_val`.
        all_explicit_imports_are_public = (;
            ignore = (
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqNewtonAlgorithm, Symbol("@cache"), :alg_cache,
                :perform_step!, :get_fsalfirstlast, :generic_solver_docstring,
                :issplit, :_fixup_ad, :_unwrap_val,
                :NLNewton, :build_nlsolver, :markfirststage!, :nlsolve!,
                :nlsolvefail, :du_alias_or_new,
            )),
    ),
)

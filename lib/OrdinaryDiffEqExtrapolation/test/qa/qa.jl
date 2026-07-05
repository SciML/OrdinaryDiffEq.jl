using SciMLTesting, OrdinaryDiffEqExtrapolation, Test

# After the solver-author extension API was declared `public` in OrdinaryDiffEqCore /
# OrdinaryDiffEqDifferentiation / OrdinaryDiffEqNonlinearSolve / DiffEqBase, the only
# remaining ExplicitImports public-API exceptions are the genuinely-internal names
# below that have no public replacement: OrdinaryDiffEqCore's private `@threaded`
# codegen macro and `_resolved_QT` controller helper, plus a handful of upstream
# (SciMLBase / FastPower / Base.Threads) internals. Everything else is clean.
# Tracked for a future make-public pass; see SciML/OrdinaryDiffEq.jl#3776.
run_qa(
    OrdinaryDiffEqExtrapolation;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore private codegen macro (deliberately kept non-public)
                Symbol("@threaded"),
                # SciMLBase internals (reshaping / val-unwrap helpers, no public replacement)
                :_reshape, :_unwrap_val, :_vec,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal (controller QT resolution)
                :_resolved_QT,
                # other upstream internals
                :fastpower,     # FastPower
                :has_Wfact,     # SciMLBase
                :maxthreadid,   # Base.Threads
            ),
        ),
    ),
)

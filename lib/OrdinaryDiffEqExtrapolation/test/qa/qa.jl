using SciMLTesting, OrdinaryDiffEqExtrapolation, Test

# After the solver-author extension API was declared `public` in OrdinaryDiffEqCore /
# OrdinaryDiffEqDifferentiation / DiffEqBase, the only remaining ExplicitImports
# public-API exceptions are the genuinely-internal names below that have no public
# replacement yet. They are imported/accessed from their true owner module, so only
# the two public-API checks flag them; everything else is clean.
# Each name is tracked for a future make-public pass; see SciML/OrdinaryDiffEq.jl#3776.
run_qa(
    OrdinaryDiffEqExtrapolation;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal (codegen macro + perf/docstring/order helpers)
                Symbol("@threaded"), :CompiledFloats, :_fixup_ad,
                :differentiation_rk_docstring, :get_current_alg_order,
                :reset_alg_dependent_opts!,
                # SciMLBase internals (reshaping / val-unwrap helpers)
                :_reshape, :_unwrap_val, :_vec,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal (controller QT resolution)
                :_resolved_QT, :resolve_basic,
                # other upstream internals
                :fastpower,     # FastPower
                :has_Wfact,     # SciMLBase
                :maxthreadid,   # Base.Threads
            ),
        ),
    ),
)

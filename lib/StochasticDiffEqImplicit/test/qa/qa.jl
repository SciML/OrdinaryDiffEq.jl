using SciMLTesting, StochasticDiffEqImplicit, Test
using JET

run_qa(
    StochasticDiffEqImplicit;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` is owned by FastBroadcast and re-exported through DiffEqBase; FastBroadcast
        # is not a direct dependency, so it is imported from DiffEqBase by design.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        # Internal (non-`public`) names from SciML packages used to implement the solvers.
        # These have no public replacement yet; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # `@..` (FastBroadcast macro re-exported via DiffEqBase) is not `public` there.
                Symbol("@.."),
                # StochasticDiffEqCore internals
                Symbol("@cache"),
                # OrdinaryDiffEqCore internals
                :perform_step!, :issplit, :default_controller, :PIController,
                :current_extrapolant!, :isnewton, :set_new_W!, :get_W, :_fixup_ad,
                # OrdinaryDiffEqNonlinearSolve internals
                :nlsolvefail, :nlsolve!, :build_nlsolver, :markfirststage!, :NLNewton,
                # OrdinaryDiffEqDifferentiation internals
                :calc_J, :calc_J!, :dolinsolve,
                # DiffEqBase internals
                :calculate_residuals, :calculate_residuals!,
                # SciMLBase internals (owner-correct but not yet `public`)
                :has_Wfact, :_vec, :_reshape, :_unwrap_val,
            )),
        # SciMLBase / OrdinaryDiffEqCore internals accessed via qualified name.
        all_qualified_accesses_are_public = (;
            ignore = (:AlgorithmInterpretation, :alg_interpretation, :set_EEst!)),
    ),
)

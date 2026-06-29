using SciMLTesting, StochasticDiffEqLeaping, Test
using JET

run_qa(
    StochasticDiffEqLeaping;
    jet_kwargs = (; target_defined_modules = true),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` is the SciML fused-broadcast macro, owned by FastBroadcast and
        # re-exported through DiffEqBase (FastBroadcast is not a direct dep here).
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        # `pois_rand` is owned by PoissonRandom and re-exported by JumpProcesses,
        # the natural dependency that owns the jump-process interface here.
        all_qualified_accesses_via_owners = (; ignore = (:pois_rand,)),
        # Internal (non-public) solver-interface names this method package must
        # extend/import. These are private API of OrdinaryDiffEqCore,
        # OrdinaryDiffEqNonlinearSolve, StochasticDiffEqCore, and DiffEqBase; they
        # have no public replacement (track make-public in SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@.."), Symbol("@cache"), :DummyController, :NLFunctional,
                :issplit, :perform_step!, :step_accept_controller!,
                :step_reject_controller!, :stepsize_controller!,
            )
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :ODEIntegrator, :build_nlsolver, :default_controller, :get_EEst,
                :isaposteriori, :markfirststage!, :nlsolve!, :nlsolve_f,
                :nlsolvefail, :pois_rand, :set_EEst!,
            )
        ),
    ),
)

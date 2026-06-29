using SciMLTesting, ImplicitDiscreteSolve, Test

run_qa(
    ImplicitDiscreteSolve;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        # Non-public OrdinaryDiffEqCore internals we extend/dispatch on. These are
        # the solver-interface generics and abstract types that OrdinaryDiffEqCore
        # has not (yet) marked `public`; tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqMutableCache, :alg_cache,
                :get_fsalfirstlast, :isfsal, :perform_step!, :isdiscretecache,
                :isdiscretealg, :beta2_default, :beta1_default, :dt_required,
                :_initialize_dae!, :allows_null_u0, :AbstractController,
                :AbstractControllerCache,
            ),
        ),
        # Non-public names accessed by qualification: the OrdinaryDiffEqCore
        # controller interface, NonlinearSolveBase solver internals, and
        # SciMLBase.has_initializeprob. Same make-public tracking as above.
        all_qualified_accesses_are_public = (;
            ignore = (
                :CommonControllerOptions, :_resolved_QT, :accept_step_controller,
                :default_controller, :resolve_basic, :setup_controller_cache,
                :step_accept_controller!, :step_reject_controller!,
                :stepsize_controller!, :sync_controllers!,
                :get_fu, :get_u, :not_terminated, :update_from_termination_cache!,
                :update_trace!, :has_initializeprob,
            ),
        ),
    ),
)

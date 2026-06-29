using SciMLTesting, OrdinaryDiffEqNordsieck, Test

run_qa(
    OrdinaryDiffEqNordsieck;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` is the package's public-API surface; the
        # reexport idiom is an unavoidable implicit import of these three names.
        no_implicit_imports = (;
            ignore = (:SciMLBase, :Reexport, Symbol("@reexport")),
        ),
        # A Nordsieck solver sublibrary necessarily extends and consumes the
        # OrdinaryDiffEqCore solver interface (cache/algorithm abstract types,
        # `alg_cache`/`perform_step!`/controller hooks, error-estimate accessors,
        # docstring helper, interpolation primitives), plus the DiffEqBase
        # `calculate_residuals`/`calculate_residuals!` residual primitives and the
        # OrdinaryDiffEqTsit5 first-step cache types. None of these are declared
        # `public`/exported by their owners yet; drop each entry as it is made
        # public upstream (tracked in SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal solver interface
                :AbstractController, :AbstractControllerCache,
                :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :alg_adaptive_order, :alg_cache, :perform_step!,
                :stepsize_controller!, :step_accept_controller!,
                :step_reject_controller!, :accept_step_controller,
                :post_newton_controller!, :default_controller,
                :setup_controller_cache, :CommonControllerOptions,
                :resolve_basic, :_resolved_QT,
                :get_current_alg_order, :get_current_adaptive_order,
                :get_fsalfirstlast, :get_EEst, :set_EEst!, :increment_nf!,
                :get_qmin, :get_qmax, :get_qsteady_min, :get_qsteady_max,
                :qmin_default, :qmax_default,
                :qsteady_min_default, :qsteady_max_default,
                :generic_solver_docstring, :trivial_limiter!,
                :ode_interpolant, :ode_interpolant!,
                # DiffEqBase residual primitives
                :calculate_residuals, :calculate_residuals!,
                # OrdinaryDiffEqTsit5 first-step caches
                :Tsit5Cache, :Tsit5ConstantCache,
            ),
        ),
    ),
)

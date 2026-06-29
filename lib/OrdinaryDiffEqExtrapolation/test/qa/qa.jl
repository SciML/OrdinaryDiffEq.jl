using SciMLTesting, OrdinaryDiffEqExtrapolation, Test

# The only remaining ExplicitImports exceptions are genuinely-internal (non-`public`)
# names of upstream SciML packages that have no public replacement yet. They are
# imported/accessed from their true owner module, so only the two public-API checks
# still flag them; everything else (implicit imports, stale imports, owner) is clean.
# Each name below is tracked for a future make-public pass; see SciML/OrdinaryDiffEq.jl#3776.
run_qa(
    OrdinaryDiffEqExtrapolation;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals (solver/cache/controller interface)
                Symbol("@cache"), Symbol("@threaded"),
                :AbstractController, :AbstractControllerCache, :CompiledFloats,
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :PIController, :PolyesterThreads, :_fixup_ad, :accept_step_controller,
                :alg_cache, :alg_maximum_order, :beta1_default, :beta2_default,
                :constvalue, :default_controller, :differentiation_rk_docstring,
                :gamma_default, :generic_solver_docstring, :get_current_adaptive_order,
                :get_current_alg_order, :get_fsalfirstlast, :get_qmax, :get_qmin,
                :isthreaded, :perform_step!, :qmin_default, :reset_alg_dependent_opts!,
                :setup_controller_cache, :step_accept_controller!, :step_reject_controller!,
                :stepsize_controller!, :unwrap_alg,
                # OrdinaryDiffEqDifferentiation internals (Jacobian/W construction)
                :build_grad_config, :build_jac_config, :calc_J, :calc_J!,
                :dolinsolve, :jacobian2W!,
                # SciMLBase internals (derivative/jacobian wrappers, reshaping helpers)
                :TimeDerivativeWrapper, :TimeGradientWrapper, :UDerivativeWrapper,
                :UJacobianWrapper, :_reshape, :_unwrap_val, :_vec,
                # DiffEqBase internals (residual/dt-floor helpers)
                :calculate_residuals, :calculate_residuals!, :timedepentdtmin,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals accessed via qualified `OrdinaryDiffEqCore.x`
                :CommonControllerOptions, :_resolved_QT, :get_EEst, :increment_nf!,
                :resolve_basic, :set_EEst!, :sync_controllers!,
                # other upstream internals
                :fastpower,     # FastPower
                :has_Wfact,     # SciMLBase
                :maxthreadid,   # Base.Threads
            ),
        ),
    ),
)

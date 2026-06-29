using SciMLTesting, OrdinaryDiffEqFIRK, Test

run_qa(
    OrdinaryDiffEqFIRK;
    explicit_imports = true,
    ei_kwargs = (
        # These names are accessed/imported from their true owner module but are
        # not (yet) declared `public`/exported there. They are internal solver-framework
        # entry points of the OrdinaryDiffEq stack (and a few helper deps) that FIRK
        # legitimately extends/uses; nothing here is fixable from this sublibrary.
        # See SciML/OrdinaryDiffEq.jl#3776 for the make-public tracking.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal integrator-stats / error-estimate API
                :get_EEst, :set_EEst!, :increment_nf!,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal solver-framework API
                :unwrap_alg, :default_controller, :PredictiveController, :PIController,
                :OrdinaryDiffEqNewtonAdaptiveAlgorithm, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqConstantCache, :alg_cache, Symbol("@threaded"),
                :isthreaded, :constvalue, :differentiation_rk_docstring,
                :trivial_limiter!, :qmax_default, :alg_adaptive_order,
                :step_accept_controller!, :step_reject_controller!, :alg_can_repeat_jac,
                :NewtonAlgorithm, :get_current_adaptive_order, :get_fsalfirstlast,
                :get_current_alg_order, :isfirk, :_fixup_ad, :perform_step!,
                :LinearAliasSpecifier, :set_discontinuity, :get_current_qmax,
                :PredictiveControllerCache, :_ode_interpolant, :_ode_interpolant!,
                :_ode_addsteps!, :has_stiff_interpolation,
                # OrdinaryDiffEqCore-owned convergence states / cutoff (re-exported by NonlinearSolve)
                :Convergence, :FastConvergence, :NLStatus, :VerySlowConvergence,
                :Divergence, :get_new_W_γdt_cutoff,
                # SciMLBase internal helpers / wrappers
                :_vec, :_reshape, :_unwrap_val, :UDerivativeWrapper, :UJacobianWrapper,
                # DiffEqBase internal residual / value helpers
                :calculate_residuals, :calculate_residuals!, :value,
                # OrdinaryDiffEqDifferentiation internal Jacobian / linsolve API
                :build_J_W, :build_jac_config, :calc_J!, :calc_J, :dolinsolve,
                :islinearfunction,
                # FastPower / SciMLOperators internal names
                :fastpower, :AbstractSciMLOperator,
            ),
        ),
    ),
)

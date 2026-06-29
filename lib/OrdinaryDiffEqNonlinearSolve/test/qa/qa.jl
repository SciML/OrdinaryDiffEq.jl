using SciMLTesting, OrdinaryDiffEqNonlinearSolve, Test

run_qa(
    OrdinaryDiffEqNonlinearSolve;
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported from a re-exporter rather than their defining package.
        # `@SciMLMessage` is owned by SciMLLogging but reached through OrdinaryDiffEqCore
        # (SciMLLogging is not a direct dependency); `WOperator`/`StaticWOperator` are the
        # OrdinaryDiffEqDifferentiation W-operator types, which ExplicitImports attributes
        # to SciMLOperators (their abstract supertype's package).
        all_explicit_imports_via_owners = (;
            ignore = (Symbol("@SciMLMessage"), :WOperator, :StaticWOperator),
        ),
        # Non-public names of upstream packages with no public alternative. Each is owned
        # by the package it is imported from; remove an entry once that package marks the
        # name `public`/exports it. Tracked in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals (solver interface)
                :AbstractNLSolver, :AbstractNLSolverAlgorithm, :AbstractNLSolverCache,
                :NewtonAlgorithm, :DAEAlgorithm, :TryAgain, :DIRK, :COEFFICIENT_MULTISTEP,
                :Convergence, :Divergence, :NLStatus, :MethodType, :error_constant,
                :alg_extrapolates, :resize_J_W!, :alg_autodiff, :has_special_newton_error,
                :_initialize_dae!, :resize_nlsolver!, :nlsolve_f, :set_new_W!, :set_W_γdt!,
                :default_nlsolve, :isnewton, :get_W, :isfirstcall, :isfirststage,
                :isJcurrent, :get_new_W_γdt_cutoff, :apply_step!, Symbol("@SciMLMessage"),
                # SciMLBase internals (owner of these names but not public)
                :_vec, :_reshape, :postamble!,
                # DiffEqBase internals
                :OrdinaryDiffEqTag, :calculate_residuals, :calculate_residuals!,
                # SciMLOperators abstract type (not public)
                :AbstractSciMLOperator,
                # OrdinaryDiffEqDifferentiation internals (Jacobian/W machinery)
                :update_W!, :is_always_new, :build_uf, :build_J_W, :WOperator,
                :StaticWOperator, :wrapprecs, :build_jac_config, :dolinsolve,
                :resize_jac_config!, :jacobian2W!, :jacobian!,
                # ForwardDiff / StaticArraysCore internals
                :Dual, :StaticArray,
            ),
        ),
        # Qualified `Module.name` accesses to non-public names with no public alternative.
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase internals (`postamble!` is extended via `SciMLBase.postamble!`)
                :value, :anyeltypedual, :FullSpecialize, :forwarddiff_chunksize,
                :has_Wfact, :has_Wfact_t, :has_jac_u, :has_jac_du, :postamble!,
                # OrdinaryDiffEqCore internals
                :ODEIntegrator, :increment_nf!,
                # ForwardDiff internals
                :Tag, :pickchunksize,
            ),
        ),
    ),
)

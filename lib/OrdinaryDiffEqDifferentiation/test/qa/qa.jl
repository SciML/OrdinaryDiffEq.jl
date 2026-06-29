using SciMLTesting, OrdinaryDiffEqDifferentiation, Test

# Residual ExplicitImports ignores: every name below is a genuinely non-public
# (or non-owner re-exported) symbol of another package that this differentiation
# sublibrary legitimately relies on. The OrdinaryDiffEqCore / SciMLBase / DiffEqBase
# entries are the internal solver/W-matrix/Jacobian-interface API that the sublib is
# built against and are make-public candidates upstream (tracked in SciML/OrdinaryDiffEq.jl#3776).
run_qa(
    OrdinaryDiffEqDifferentiation;
    aqua_kwargs = (; piracies = false, ambiguities = false),
    explicit_imports = true,
    ei_kwargs = (;
        # `@set` (Accessors) and `@SciMLMessage` (SciMLLogging) are re-exported by
        # SciMLBase / OrdinaryDiffEqCore; their owners are not direct deps.
        all_explicit_imports_via_owners = (; ignore = (Symbol("@set"), Symbol("@SciMLMessage"))),
        all_qualified_accesses_via_owners = (; ignore = (Symbol("@set"),)),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal solver API (make-public candidates)
                :AbstractNLSolver, :alg_autodiff, :CompositeAlgorithm, :concrete_jac,
                :constvalue, :DAEAlgorithm, :diffdir, :Divergence, :get_chunksize,
                :get_new_W_γdt_cutoff, :get_W, :isfirstcall, :isfirststage, :isJcurrent,
                :isnewton, :issplit, :isWmethod, :nlsolve_f,
                :OrdinaryDiffEqAdaptiveExponentialAlgorithm,
                :OrdinaryDiffEqAdaptiveImplicitAlgorithm, :OrdinaryDiffEqAlgorithm,
                :OrdinaryDiffEqCache, :OrdinaryDiffEqExponentialAlgorithm,
                :OrdinaryDiffEqImplicitAlgorithm, :resize_J_W!, :set_new_W!,
                Symbol("set_W_γdt!"), :TryAgain, :unwrap_alg,
                # SciMLBase internal wrappers / helpers (make-public candidates)
                :UDerivativeWrapper, :UJacobianWrapper, :_unwrap_val, :_vec,
                # other framework internals
                :AbstractSciMLOperator,  # SciMLOperators
                :OrdinaryDiffEqTag,      # DiffEqBase
                :StaticArray, :StaticMatrix,  # StaticArraysCore
                # Accessors `@set` / SciMLLogging `@SciMLMessage` re-exported by
                # SciMLBase / OrdinaryDiffEqCore (owners are not direct deps)
                Symbol("@set"), Symbol("@SciMLMessage"),
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base / framework internals accessed by qualification
                Symbol("@pure"),         # Base
                Symbol("@set"),          # Accessors macro re-exported by SciMLBase
                :AbstractSciMLOperator,  # SciMLOperators
                # LinearSolve internals
                :DefaultLinearSolver, :InvPreconditioner, :LinearCache,
                :init_cacheval, :needs_concrete_A,
                # ForwardDiff internals
                :JacobianConfig, :Tag, :pickchunksize,
                # OrdinaryDiffEqCore internal solver API (make-public candidates)
                :_get_fwd_chunksize_int, :get_EEst, Symbol("increment_nf!"), :unwrap_alg,
                # DiffEqBase internals
                :default_factorize, :prepare_alg,
                # SciMLBase Jacobian-interface predicates (make-public candidates)
                :has_Wfact_t, :has_colorvec, :has_jac_du, :has_jac_u,
                # DifferentiationInterface internals
                Symbol("prepare!_derivative"), Symbol("prepare!_jacobian"),
                # this sublib's own internal sparse-handling API, accessed
                # qualified from OrdinaryDiffEqDifferentiationSparseArraysExt
                :get_nzval, :is_sparse, :is_sparse_csc, :nonzeros,
                Symbol("set_all_nzval!"), :spzeros,
            ),
        ),
    ),
)

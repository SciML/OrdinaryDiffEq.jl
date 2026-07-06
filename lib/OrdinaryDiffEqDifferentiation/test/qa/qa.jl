using SciMLTesting, OrdinaryDiffEqDifferentiation, Test

# Residual ExplicitImports ignores. After the solver-author extension API of
# OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation / DiffEqBase was declared
# `public` (Julia 1.11+ `public`), the only names left are genuine non-public
# externals, this sublib's own internal sparse-handling API, and a small set of
# OrdinaryDiffEqCore names not yet in that public block (make-public follow-ups,
# tracked in SciML/OrdinaryDiffEq.jl#3776).
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
                # Accessors `@set` / SciMLLogging `@SciMLMessage` re-exported by
                # SciMLBase / OrdinaryDiffEqCore (owners are not direct deps)
                Symbol("@set"), Symbol("@SciMLMessage"),
                # SciMLBase internal helpers (make-public candidates)
                :_unwrap_val, :_vec,
                # OrdinaryDiffEqCore internals not yet in its public block
                :concrete_jac, :diffdir, :get_chunksize, :isnewton,
                # other framework internals
                :AbstractSciMLOperator,       # SciMLOperators
                :StaticArray, :StaticMatrix,  # StaticArraysCore
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                Symbol("@pure"),         # Base
                Symbol("@set"),          # Accessors macro re-exported by SciMLBase
                :AbstractSciMLOperator,  # SciMLOperators
                # LinearSolve internals
                :DefaultLinearSolver, :InvPreconditioner, :LinearCache,
                :init_cacheval, :needs_concrete_A,
                # ForwardDiff internals
                :JacobianConfig, :Tag, :pickchunksize,
                # OrdinaryDiffEqCore internal not yet in its public block
                :_get_fwd_chunksize_int,
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

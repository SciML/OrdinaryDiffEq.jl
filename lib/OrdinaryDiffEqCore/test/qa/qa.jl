using SciMLTesting, OrdinaryDiffEqCore, Test

# `OrdinaryDiffEqCore.Predictor` is a generated `EnumX.@enumx` submodule, which
# ExplicitImports cannot statically analyze; allow it to be unanalyzable.
const UNANALYZABLE = (OrdinaryDiffEqCore.Predictor,)

run_qa(
    OrdinaryDiffEqCore;
    aqua_kwargs = (; piracies = false, unbound_args = false),
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (; allow_unanalyzable = UNANALYZABLE),
        # These names are not used by OrdinaryDiffEqCore itself, but are imported
        # into its namespace and re-imported by dependent OrdinaryDiffEq.jl
        # sublibraries (e.g. `import OrdinaryDiffEqCore: calculate_residuals!`).
        # ExplicitImports cannot see that cross-package usage, so it reports them
        # as stale; they must remain part of this package's namespace contract.
        no_stale_explicit_imports = (;
            allow_unanalyzable = UNANALYZABLE,
            ignore = (
                :BrownFullBasicInit, :ShampineCollocationInit, :DEVerbosity,
                :Minimal, :_vec, :_reshape, :unwrap_cache,
                :calculate_residuals, :calculate_residuals!,
            ),
        ),
        # `constructorof` is owned by ConstructionBase but reached through SciMLBase
        # (not a direct dependency); accessing it via SciMLBase is intentional.
        all_qualified_accesses_via_owners = (; ignore = (:constructorof,)),
        # Internal (non-`public`) names of upstream packages that OrdinaryDiffEqCore
        # genuinely needs and that have no public replacement yet.
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base / Core internals
                Symbol("@max_methods"), :Experimental, :Typeof, :promote_op,
                # ConstructionBase (owned there, accessed via SciMLBase)
                :constructorof,
                # SciMLOperators internal
                :AbstractSciMLOperator,
                # EnzymeCore / EnzymeCore.EnzymeRules internals
                :EnzymeRules, :inactive_noinl,
                # SciMLBase internals with no public replacement yet
                :enable_interpolation_sensitivitymode,
                :forwarddiff_chunksize, :get_root_indp,
                :get_save_idxs_and_saved_subsystem, :has_initializeprob,
                :has_lazy_interpolation, :late_binding_update_u0_p, :remaker_of,
                :save_discretes_if_enabled!, :save_final_discretes!,
                :strip_interpolation, :struct_as_namedtuple, :unitfulvalue,
                :unwrap_cache, :value,
            ),
        ),
        # Internal (non-`public`) names imported from upstream packages.
        all_explicit_imports_are_public = (;
            ignore = (
                # TruncatedStacktraces internals
                Symbol("@truncate_stacktrace"), :VERBOSE_MSG,
                # FunctionWrappers internal
                :FunctionWrapper,
                # FastPower internal
                :fastpower,
                # SciMLBase internals
                :SENSITIVITY_INTERP_MESSAGE, :_unwrap_val, :last_step_failed,
                :postamble!,
                # DiffEqBase internal
                :_process_verbose_param,
                # Non-public names re-exported for dependent sublibraries
                # (see `no_stale_explicit_imports` above).
                :_vec, :_reshape, :unwrap_cache,
            ),
        ),
    ),
)

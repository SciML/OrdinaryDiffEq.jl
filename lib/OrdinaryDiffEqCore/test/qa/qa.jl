using SciMLTesting, OrdinaryDiffEqCore, Test

# `OrdinaryDiffEqCore.Predictor` is a generated `EnumX.@enumx` submodule, which
# ExplicitImports cannot statically analyze; allow it to be unanalyzable.
const UNANALYZABLE = (OrdinaryDiffEqCore.Predictor,)

@static if VERSION >= v"1.11.0-DEV.469"
    for name in (
            :AbstractNLSolverCache, :AutoSwitchCache, :DefaultCache,
            :handle_callback_modifiers!, :isdiscretecache,
            :resolve_stage_step_limiters, :trivial_limiter!,
        )
        @test Base.ispublic(OrdinaryDiffEqCore, name)
    end
end

run_qa(
    OrdinaryDiffEqCore;
    aqua_kwargs = (; piracies = false, unbound_args = false),
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
                :Minimal, :_vec, :_reshape, :unwrap_cache, :_unwrap_val,
                :calculate_residuals, :calculate_residuals!,
            ),
        ),
        # Internal (non-`public`) names of upstream packages that OrdinaryDiffEqCore
        # genuinely needs and that have no public replacement yet.
        all_qualified_accesses_are_public = (;
            ignore = (
                # DiffEqBase — owner-internal, no public alternative
                :NAN_CHECK,
                # Base / Core internals
                Symbol("@max_methods"), :Experimental, :Typeof, :promote_op,
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
                # DiffEqBase — owner-internal, no public alternative
                :NAN_CHECK,
                # TruncatedStacktraces internals
                Symbol("@truncate_stacktrace"), :VERBOSE_MSG,
                # FunctionWrappers internal
                :FunctionWrapper,
                # FastPower internal
                :fastpower,
                # SciMLBase internals
                :SENSITIVITY_INTERP_MESSAGE, :last_step_failed,
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

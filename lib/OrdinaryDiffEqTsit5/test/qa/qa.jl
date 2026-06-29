using SciMLTesting, OrdinaryDiffEqTsit5, Test

run_qa(
    OrdinaryDiffEqTsit5;
    explicit_imports = true,
    ei_kwargs = (;
        # Solver-internal names accessed via their (qualified) owner module that
        # have no public alternative on the registered releases. The OrdinaryDiffEqCore
        # step-statistics mutators (increment_nf!/set_EEst!) and the lorenz/lorenz_oop
        # precompile-workload test functions are OrdinaryDiffEqCore internals;
        # FunctionWrapperSpecialize/NoSpecialize are SciMLBase specialization markers
        # used in the (preference-gated) precompile workload and are not declared public.
        all_qualified_accesses_are_public = (;
            ignore = (
                :increment_nf!,              # OrdinaryDiffEqCore
                :set_EEst!,                  # OrdinaryDiffEqCore
                :lorenz,                     # OrdinaryDiffEqCore
                :lorenz_oop,                 # OrdinaryDiffEqCore
                :FunctionWrapperSpecialize,  # SciMLBase
                :NoSpecialize,               # SciMLBase
            ),
        ),
        # Core solver-interface names imported from their owner module but not yet
        # declared `public`. These are the OrdinaryDiffEq internal extension/dispatch
        # API (OrdinaryDiffEqCore), the DiffEqBase residual helpers, SciMLBase's @def,
        # and the TruncatedStacktraces macro. No public alternative exists; tracked for
        # make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore
                :alg_cache, :alg_stability_size, :AutoAlgSwitch,
                Symbol("@cache"), :CompiledFloats, :CompositeAlgorithm,
                :constvalue, :DerivativeOrderNotPossibleError,
                :explicit_rk_docstring, Symbol("@fold"), :get_fsalfirstlast,
                :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!,
                Symbol("@OnDemandTableauExtract"),
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache, :perform_step!, :trivial_limiter!,
                # DiffEqBase
                :calculate_residuals, :calculate_residuals!,
                # SciMLBase
                Symbol("@def"),
                # TruncatedStacktraces
                Symbol("@truncate_stacktrace"),
            ),
        ),
    ),
)

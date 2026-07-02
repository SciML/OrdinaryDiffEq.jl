using SciMLTesting, OrdinaryDiffEqTsit5, Test

run_qa(
    OrdinaryDiffEqTsit5;
    explicit_imports = true,
    ei_kwargs = (;
        # Residual after the PHASE-A make-public sweep: every solver-author API name
        # (alg_cache/perform_step!/_ode_interpolant/... in OrdinaryDiffEqCore,
        # calculate_residuals in DiffEqBase, etc.) is now declared `public` on this
        # branch and no longer needs an ignore. What remains are owner-internal names
        # with no public alternative.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — private codegen macros kept non-public in PHASE A
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"),
                # OrdinaryDiffEqCore — precompile-workload test functions (internal)
                :lorenz, :lorenz_oop,
                # OrdinaryDiffEqCore — genuine internals (no public alternative)
                :CompiledFloats, :DerivativeOrderNotPossibleError,
                :_ode_interpolant!, :trivial_limiter!,
                # SciMLBase
                Symbol("@def"),
                # TruncatedStacktraces
                Symbol("@truncate_stacktrace"),
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — private codegen macros kept non-public in PHASE A
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"),
                # OrdinaryDiffEqCore — genuine internals (no public alternative)
                :CompiledFloats, :DerivativeOrderNotPossibleError,
                :_ode_interpolant!, :trivial_limiter!,
                # SciMLBase
                Symbol("@def"),
                # TruncatedStacktraces
                Symbol("@truncate_stacktrace"),
            ),
        ),
    ),
)

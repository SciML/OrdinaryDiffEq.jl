using SciMLTesting, OrdinaryDiffEqTsit5, Test

run_qa(
    OrdinaryDiffEqTsit5;
    explicit_imports = true,
    ei_kwargs = (;
        # Residual after the PHASE-A make-public sweep: the solver-author API
        # (alg_cache/perform_step!/_ode_interpolant!/CompiledFloats/... in
        # OrdinaryDiffEqCore, calculate_residuals in DiffEqBase, etc.) is now
        # declared `public` on this branch and no longer needs an ignore. What
        # remains are owner-internal names with no public alternative.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — precompile-workload test functions (internal)
                :lorenz, :lorenz_oop, :lorenz_p, :lorenz_p_params,
            ),
        ),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore — private codegen macros / limiter kept
                # non-public in PHASE A (owner-internal, no public alternative)
                Symbol("@fold"), Symbol("@OnDemandTableauExtract"),
                :trivial_limiter!,
                # SciMLBase
                Symbol("@def"),
                # TruncatedStacktraces
                Symbol("@truncate_stacktrace"),
            ),
        ),
    ),
)

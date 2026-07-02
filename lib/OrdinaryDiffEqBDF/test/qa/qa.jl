using SciMLTesting, OrdinaryDiffEqBDF, Test

run_qa(
    OrdinaryDiffEqBDF;
    explicit_imports = true,
    ei_kwargs = (
        # Precompile-workload test problems that live in OrdinaryDiffEqCore but are
        # deliberately kept non-public (they are test fixtures, not solver-author API).
        all_qualified_accesses_are_public = (;
            ignore = (:lorenz, :lorenz_oop),
        ),
        # Genuine solver-author internals still kept non-public by their owners.
        # Everything OrdinaryDiffEqCore made public for the solver extension API has
        # been dropped; drop each remaining entry once its owner marks it public
        # upstream (tracked in SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals (not in the public extension API)
                :DerivativeOrderNotPossibleError, :error_constant,
                :get_current_alg_order, :get_current_qmax,
                :has_special_newton_error, :has_stiff_interpolation,
                :isnewton, :qsteady_min_default, :qsteady_max_default,
                :resolve_basic, :_resolved_QT, :set_discontinuity,
                :trivial_limiter!, :_ode_interpolant!, :_fixup_ad,
                # OrdinaryDiffEqSDIRK ESDIRK-IMEX caches/tableau (owner-internal)
                :ESDIRKIMEXCache, :ESDIRKIMEXConstantCache,
                :ImplicitEulerESDIRKIMEXTableau,
                # SciMLBase internal
                :_unwrap_val,
                # TruncatedStacktraces internal macro
                Symbol("@truncate_stacktrace"),
            ),
        ),
    ),
)

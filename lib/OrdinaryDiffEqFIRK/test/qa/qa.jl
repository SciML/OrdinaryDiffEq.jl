using SciMLTesting, OrdinaryDiffEqFIRK, Test

run_qa(
    OrdinaryDiffEqFIRK;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal names deliberately kept non-public
                # (private codegen/perf macro + solver-framework internals not part of
                # the declared extension surface).
                Symbol("@threaded"), :_fixup_ad, :_ode_interpolant!,
                :alg_can_repeat_jac, :differentiation_rk_docstring,
                :get_current_alg_order, :get_current_qmax,
                :has_stiff_interpolation, :isfirk, :set_discontinuity,
                :trivial_limiter!,
                # SciMLBase internal helpers / wrappers (pending SciMLBase#1412)
                :_vec, :_reshape, :_unwrap_val, :value,
                # Genuine external deps
                :fastpower,               # FastPower
                :AbstractSciMLOperator,   # SciMLOperators
            ),
        ),
    ),
)

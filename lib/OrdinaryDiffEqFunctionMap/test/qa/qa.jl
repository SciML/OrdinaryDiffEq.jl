using SciMLTesting, OrdinaryDiffEqFunctionMap, Test

run_qa(
    OrdinaryDiffEqFunctionMap;
    aqua_kwargs = (; piracies = false),  # piracy is needed for default-algorithm dispatch
    explicit_imports = true,
    ei_kwargs = (;
        # Internal (non-public) names extended/dispatched from OrdinaryDiffEqCore.
        # No public API exists for these solver-internal hooks; see SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqConstantCache, :perform_step!, :unwrap_alg,
                :alg_cache, Symbol("@cache"), :_ode_addsteps!, :_ode_interpolant,
                :_ode_interpolant!, :get_fsalfirstlast, :isfsal, :beta1_default,
                :beta2_default, :dt_required, :isdiscretecache, :isdiscretealg,
            ),
        ),
        # increment_nf! is OrdinaryDiffEqCore-internal; EvalFunc is DiffEqBase-internal;
        # DISCRETE_{IN,OUT}OFPLACE_DEFAULT are SciMLBase-owned but non-public.
        all_qualified_accesses_are_public = (;
            ignore = (
                :increment_nf!, :EvalFunc,
                :DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT,
            ),
        ),
    ),
)

using SciMLTesting, OrdinaryDiffEqExtrapolation, Test

run_qa(
    OrdinaryDiffEqExtrapolation;
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore private codegen macro (deliberately kept non-public)
                Symbol("@threaded"),
                # OrdinaryDiffEqDifferentiation owner-internal cross-sublibrary hooks;
                # no public wrapper exists.
                :build_grad_config, :build_jac_config, :calc_J, :calc_J!,
                :dolinsolve, :jacobian2W!,
                # SciMLBase internals (reshaping / val-unwrap helpers, no public replacement)
                :_reshape, :_unwrap_val, :_vec,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore owner-internal (controller QT resolution)
                :_resolved_QT,
                # other upstream internals
                :fastpower,     # FastPower
                :has_Wfact,     # SciMLBase
                :maxthreadid,   # Base.Threads
            ),
        ),
    ),
)

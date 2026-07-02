using SciMLTesting, OrdinaryDiffEqRKN, Test

run_qa(
    OrdinaryDiffEqRKN;
    explicit_imports = true,
    ei_kwargs = (;
        # Genuine solver-author internals that their owners deliberately keep
        # non-public. Everything OrdinaryDiffEqCore/DiffEqBase made public for the
        # solver extension API has been dropped; drop each remaining entry once its
        # owner marks it public upstream (tracked in SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase codegen macro (owner keeps non-public)
                Symbol("@def"),
                # DiffEqBase perf macro (owner keeps non-public)
                Symbol("@tight_loop_macros"),
                # OrdinaryDiffEqCore internals not in the public extension API
                # (the non-bang `_ode_interpolant` is public; the `!`-variant is not)
                :CompiledFloats, :_ode_interpolant!,
            ),
        ),
    ),
)

using SciMLTesting, OrdinaryDiffEqSSPRK, Test

run_qa(
    OrdinaryDiffEqSSPRK;
    explicit_imports = true,
    ei_kwargs = (
        # Residual non-public names after the solver-author API was made public
        # in OrdinaryDiffEqCore / DiffEqBase. Everything still listed here is a
        # deliberate owner-internal that is NOT part of the supported public
        # surface.
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase private codegen macro (no public counterpart).
                Symbol("@def"),
                # OrdinaryDiffEqCore codegen/perf internals kept non-public on
                # purpose (see OrdinaryDiffEqCore's public block comment):
                # the interpolation write-hook, the SSP coefficient accessor,
                # and the no-op limiter.
                :_ode_interpolant!, :ssp_coefficient, :trivial_limiter!,
            ),
        ),
        # OrdinaryDiffEqCore precompile-workload test problems (not public API).
        all_qualified_accesses_are_public = (;
            ignore = (:lorenz, :lorenz_oop),
        ),
    ),
)

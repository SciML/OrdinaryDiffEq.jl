using SciMLTesting, OrdinaryDiffEqAdamsBashforthMoulton, Test

run_qa(
    OrdinaryDiffEqAdamsBashforthMoulton;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqLowOrderRK starter-cache types (owner-internal;
                # kept non-public — used only as the multistep bootstrap cache).
                :BS3Cache, :BS3ConstantCache, :RK4Cache, :RK4ConstantCache,
                # OrdinaryDiffEqCore internal limiter hook (owner-internal).
                :trivial_limiter!,
            )),
    ),
)

using SciMLTesting, OrdinaryDiffEqSymplecticRK, Test

run_qa(
    OrdinaryDiffEqSymplecticRK;
    explicit_imports = true,
    ei_kwargs = (;
        # OrdinaryDiffEqCore-internal (non-public) name the solver genuinely
        # needs; `CompiledFloats` is `const CompiledFloats = Union{Float32,Float64}`,
        # not declared public in Core. Tracked for make-public in
        # SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (; ignore = (:CompiledFloats,)),
    ),
)

using SciMLTesting, OrdinaryDiffEq, SciMLBase
using Test

@testset "PureKLU-compatible MuladdMacro floors" begin
    for package in (
            "OrdinaryDiffEqBDF",
            "OrdinaryDiffEqCore",
            "OrdinaryDiffEqExponentialRK",
            "OrdinaryDiffEqFunctionMap",
            "OrdinaryDiffEqMultirate",
            "OrdinaryDiffEqPDIRK",
            "OrdinaryDiffEqSDIRK",
        )
        project = read(joinpath(@__DIR__, "..", "..", "lib", package, "Project.toml"), String)
        floor_match = match(r"(?m)^MuladdMacro = \"([0-9]+\.[0-9]+\.[0-9]+)", project)
        @test floor_match !== nothing
        if floor_match !== nothing
            @test VersionNumber(floor_match[1]) >= v"0.2.4"
        end
    end
end

# The umbrella package's QA lane historically ran only the ExplicitImports checks
# (no Aqua/JET), so keep `aqua = false` to preserve that scope.
run_qa(
    OrdinaryDiffEq;
    aqua = false,
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (; skip = (Base, Core, SciMLBase)),
    ),
)

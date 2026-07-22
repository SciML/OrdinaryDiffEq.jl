using SciMLTesting, OrdinaryDiffEq, SciMLBase
using Test

@testset "PureKLU-compatible MuladdMacro floors" begin
    lib_dir = joinpath(@__DIR__, "..", "..", "lib")
    projects = filter(readdir(lib_dir)) do package
        project_path = joinpath(lib_dir, package, "Project.toml")
        isfile(project_path) || return false
        project = read(project_path, String)
        has_muladdmacro = occursin(
            r"(?m)^MuladdMacro = \"[0-9]+\.[0-9]+\.[0-9]+", project
        )
        uses_core = package == "OrdinaryDiffEqCore" ||
            occursin(r"(?m)^OrdinaryDiffEqCore = ", project)
        has_muladdmacro && uses_core
    end
    @test length(projects) == 39
    for package in projects
        project = read(joinpath(lib_dir, package, "Project.toml"), String)
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

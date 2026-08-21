using SciMLTesting, OrdinaryDiffEq
using ADTypes, CommonSolve, DiffEqBase, OrdinaryDiffEqBDF, OrdinaryDiffEqDefault,
    OrdinaryDiffEqRosenbrock, OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner,
    SciMLBase, SciMLLogging
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
    @test length(projects) == 41
    for package in projects
        project = read(joinpath(lib_dir, package, "Project.toml"), String)
        floor_match = match(r"(?m)^MuladdMacro = \"([0-9]+\.[0-9]+\.[0-9]+)", project)
        @test floor_match !== nothing
        if floor_match !== nothing
            @test VersionNumber(floor_match[1]) >= v"0.2.4"
        end
    end
end

const ORDINARYDIFFEQ_REEXPORTS = intersect(
    public_api_names(OrdinaryDiffEq),
    union(
        public_api_names(ADTypes),
        public_api_names(CommonSolve),
        public_api_names(DiffEqBase),
        public_api_names(OrdinaryDiffEqBDF),
        public_api_names(OrdinaryDiffEqDefault),
        public_api_names(OrdinaryDiffEqRosenbrock),
        public_api_names(OrdinaryDiffEqTsit5),
        public_api_names(OrdinaryDiffEqVerner),
        public_api_names(SciMLBase),
        public_api_names(SciMLLogging),
        # `successful_retcode` is intentionally available from the umbrella
        # package for the documented `solve` workflow, although older SciMLBase
        # releases expose it through `public` rather than `export`.
        (:SciMLBase, :SciMLLogging, :successful_retcode),
        # `public` on a Core-owned name counts as a public reexport; listed by
        # name so the check still rejects the rest of Core's surface.
        (:Sequential, :BaseThreads, :PolyesterThreads),
    ),
)

@testset "Re-exported SciMLBase API is uniform and follows the documented rule" begin
    # docs/src/api/reexports.md defines which SciMLBase names a solver sublibrary re-exports
    # (the user-facing common interface: problem/function types, solve, callbacks, integrator
    # interface, ensembles). That page is the single source of truth; this testset enforces it.
    repo_dir = joinpath(@__DIR__, "..", "..")

    documented = Symbol[]
    for line in eachline(joinpath(repo_dir, "docs", "src", "api", "reexports.md"))
        occursin(r"^[A-Za-z_][A-Za-z0-9_!]*(?:, [A-Za-z_][A-Za-z0-9_!]*)*,?$", line) || continue
        append!(documented, Symbol.(filter(!isempty, strip.(split(line, ",")))))
    end
    @test allunique(documented)
    @test length(documented) > 100

    exported = setdiff(names(SciMLBase), [:SciMLBase])
    @test isempty(setdiff(documented, exported))

    # Categories 2-5 are enumerated by hand, but category 1 (problem and function types) is
    # mechanical, so re-derive it: a concrete, non-dispatch-tag `...Problem`/`...Function`.
    # A new SciMLBase problem type fails here until the sublibraries re-export it.
    is_dispatch_tag(name) = startswith(name, "Abstract") || startswith(name, "Standard")
    function is_concrete_type(name)
        value = getglobal(SciMLBase, name)
        return value isa Type && !isabstracttype(Base.unwrap_unionall(value))
    end
    required = filter(exported) do name
        text = String(name)
        (endswith(text, "Problem") || endswith(text, "Function")) &&
            !is_dispatch_tag(text) && is_concrete_type(name)
    end
    @test isempty(setdiff(required, documented))

    function reexported_names(path)
        source = read(path, String)
        m = match(r"@reexport using SciMLBase:((?:[^\n]*\n)(?:    [^\n]*\n)*)", source)
        m === nothing && return nothing
        return Symbol.(filter(!isempty, strip.(split(replace(m[1], r"\s+" => " "), ","))))
    end

    # Sublibraries that still re-export all of SciMLBase are a superset of the rule and are
    # migrated separately; every sublibrary that names its re-exports must match the page.
    lib_dir = joinpath(repo_dir, "lib")
    selective = Dict{String, Vector{Symbol}}()
    for package in readdir(lib_dir)
        source = joinpath(lib_dir, package, "src", package * ".jl")
        isfile(source) || continue
        found = reexported_names(source)
        found === nothing || (selective[package] = found)
    end
    @test length(selective) == 9
    for (package, found) in selective
        @test (package, sort(found)) == (package, sort(documented))
    end
end

# The umbrella package's QA lane historically ran only the ExplicitImports checks
# (no Aqua/JET), so keep `aqua = false` to preserve that scope.
run_qa(
    OrdinaryDiffEq;
    aqua = false,
    explicit_imports = true,
    check_reexports = true,
    reexports_allow = ORDINARYDIFFEQ_REEXPORTS,
    ei_kwargs = (;
        no_implicit_imports = (; skip = (Base, Core, SciMLBase)),
    ),
)

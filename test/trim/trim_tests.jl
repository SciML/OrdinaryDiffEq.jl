# Trim (static-compilation) tests for OrdinaryDiffEq solver sub-packages.
#
# Each solver is compiled into a standalone native executable with
# JuliaC --trim=safe, then run and validated against known output.
#
# To add a new solver test, add a SolverConfig entry to SOLVER_CONFIGS.

using Test
using JuliaC

@assert VERSION >= v"1.12.0" "Trim tests require Julia 1.12+"

const TRIM_DIR = @__DIR__
const LOG_DIR = joinpath(TRIM_DIR, "logs")

# ── Solver configuration ──────────────────────────────────────────────

Base.@kwdef struct SolverConfig
    pkg::String           # OrdinaryDiffEq sub-package exporting the solver
    uuid::String          # that package's UUID
    alg_type::String      # the algorithm struct name to import
    constructor::String   # Julia expression constructing the solver instance
    display_name::String  # human-readable name for log messages
end

const SOLVER_CONFIGS = Dict{String,SolverConfig}(
    "tsit5" => SolverConfig(
        pkg = "OrdinaryDiffEqTsit5",
        uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a",
        alg_type = "Tsit5",
        constructor = "Tsit5()",
        display_name = "Tsit5",
    ),
    "fbdf" => SolverConfig(
        pkg = "OrdinaryDiffEqBDF",
        uuid = "6ad6398a-0878-4a85-9266-38940aa047c8",
        alg_type = "FBDF",
        constructor = "FBDF(autodiff = AutoForwardDiff(chunksize = 1), linsolve = LUFactorization())",
        display_name = "FBDF",
    ),
    "rodas5p" => SolverConfig(
        pkg = "OrdinaryDiffEqRosenbrock",
        uuid = "43230ef6-c299-4910-a778-202eb28ce4ce",
        alg_type = "Rodas5P",
        constructor = "Rodas5P(autodiff = AutoForwardDiff(chunksize = 1), linsolve = LUFactorization())",
        display_name = "Rodas5P",
    ),
    "auto" => SolverConfig(
        pkg = "OrdinaryDiffEqDefault",
        uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7",
        alg_type = "DefaultODEAlgorithm",
        constructor = "DefaultODEAlgorithm(autodiff = AutoForwardDiff(chunksize = 1), linsolve = LUFactorization())",
        display_name = "DefaultODEAlgorithm",
    ),
)

const SOLVER_ORDER = ["tsit5", "fbdf", "rodas5p", "auto"]

# Path from a generated sub-project back to the monorepo lib/ directory.
# Generated projects live at test/trim/<solver>/, so ../../.. reaches the repo root.
const LIB_REL_PATH = "../../.."

# All deps for generated sub-projects: name => (uuid, lib_subdir_or_nothing).
# Packages with a non-nothing lib_subdir get a [sources] entry pointing to
# the local monorepo checkout.
const ALL_DEPS = Dict(
    "ADTypes"                  => ("47edcb42-4c32-4615-8424-f2b9edc5f35b", nothing),
    "LinearSolve"              => ("7ed4a6bd-45f5-4d41-b270-4a48e9bafcae", nothing),
    "SciMLBase"                => ("0bca4576-84f4-4d90-8ffe-ffa030f20462", nothing),
    "SciMLLogging"             => ("a6db7da4-7206-11f0-1eab-35f2a5dbe1d1", nothing),
    "OrdinaryDiffEqCore"       => ("bbf590c4-e513-4bbe-9b18-05decba2e5d8", "OrdinaryDiffEqCore"),
    "OrdinaryDiffEqTsit5"      => ("b1df2697-797e-41e3-8120-5422d3b24e4a", "OrdinaryDiffEqTsit5"),
    "OrdinaryDiffEqBDF"        => ("6ad6398a-0878-4a85-9266-38940aa047c8", "OrdinaryDiffEqBDF"),
    "OrdinaryDiffEqRosenbrock" => ("43230ef6-c299-4910-a778-202eb28ce4ce", "OrdinaryDiffEqRosenbrock"),
    "OrdinaryDiffEqDefault"    => ("50262376-6c5a-4cf5-baba-aaf4f84d72d7", "OrdinaryDiffEqDefault"),
    "DiffEqBase"               => ("2b5f629d-d688-5b77-993f-72d75c75574e", "DiffEqBase"),
)

# ── Code generation ────────────────────────────────────────────────────

function generate_project_toml(cfg::SolverConfig)
    deps = Dict{String,String}()
    sources = Dict{String,String}()

    # Add all shared deps + local packages
    for (pkg, (uuid, lib_subdir)) in ALL_DEPS
        deps[pkg] = uuid
        if lib_subdir !== nothing
            sources[pkg] = "$LIB_REL_PATH/lib/$lib_subdir"
        end
    end
    # The solver's own package (may already be in ALL_DEPS, that's fine — Dict deduplicates)
    deps[cfg.pkg] = cfg.uuid

    deps_str = join(["$k = \"$v\"" for (k, v) in sort(collect(deps))], "\n")
    sources_str = join(["$k = {path = \"$v\"}" for (k, v) in sort(collect(sources))], "\n")

    return """
    [deps]
    $deps_str

    [sources]
    $sources_str
    """
end

function generate_main_jl(name::String, cfg::SolverConfig)
    return """
    using $(cfg.pkg): $(cfg.alg_type)
    using SciMLBase: ODEProblem, ODEFunction, solve
    using ADTypes: AutoForwardDiff
    using LinearSolve: LUFactorization
    using OrdinaryDiffEqCore: DEVerbosity
    using SciMLLogging: None

    function lotka_volterra!(du, u, p, t)::Nothing
        α, β, γ, δ = p
        du[1] = α * u[1] - β * u[1] * u[2]
        du[2] = δ * u[1] * u[2] - γ * u[2]
        return nothing
    end

    function run_solve(outfile::String, α::Float64, β::Float64)::Int
        p = [α, β, 2.0, 1.0]
        u0 = [1.0, 1.0]
        tspan = (0.0, 10.0)

        f = ODEFunction{true}(lotka_volterra!)
        prob = ODEProblem(f, u0, tspan, p)
        solver = $(cfg.constructor)
        sol = solve(prob, solver; saveat = 0.1, abstol = 1e-8, reltol = 1e-8,
                    verbose = DEVerbosity(None()))

        io = open(outfile, "w")
        for i in eachindex(sol.t)
            print(io, sol.t[i], ",", sol.u[i][1], ",", sol.u[i][2], "\\n")
        end
        close(io)

        return length(sol.t)
    end

    function (@main)(args::Vector{String})::Int32
        if length(args) < 1
            Core.print(Core.stderr, "Usage: $(name) <outfile> [α] [β]\\n")
            return Int32(1)
        end
        outfile = args[1]
        α = length(args) >= 2 ? parse(Float64, args[2]) : 1.5
        β = length(args) >= 3 ? parse(Float64, args[3]) : 1.0

        n = run_solve(outfile, α, β)
        Core.println("Solved Lotka-Volterra with $(cfg.display_name): \$(n) timesteps written to \$(outfile)")
        return Int32(0)
    end
    """
end

function generate_test_project(name::String, cfg::SolverConfig)
    test_dir = joinpath(TRIM_DIR, name)
    mkpath(test_dir)
    write(joinpath(test_dir, "Project.toml"), generate_project_toml(cfg))
    write(joinpath(test_dir, "main.jl"), generate_main_jl(name, cfg))
    # Remove stale Manifest so Pkg.instantiate picks up workspace changes
    manifest = joinpath(test_dir, "Manifest.toml")
    isfile(manifest) && rm(manifest)
    return test_dir
end

# ── Compile & validate ─────────────────────────────────────────────────

function compile_solver(solver_name::String, build_dir::String)
    test_dir = joinpath(TRIM_DIR, solver_name)
    main_jl = joinpath(test_dir, "main.jl")

    outname = joinpath(build_dir, solver_name)
    logfile = joinpath(LOG_DIR, "$(solver_name)_compile.log")

    open(logfile, "w") do log
        redirect_stdio(; stdout = log, stderr = log) do
            img = JuliaC.ImageRecipe(;
                output_type = "--output-exe",
                trim_mode = "safe",
                file = abspath(main_jl),
                project = abspath(test_dir),
                verbose = true,
            )
            JuliaC.compile_products(img)

            link = JuliaC.LinkRecipe(; image_recipe = img, outname = outname)
            JuliaC.link_products(link)

            bundle = JuliaC.BundleRecipe(; link_recipe = link, output_dir = build_dir)
            JuliaC.bundle_products(bundle)
        end
    end

    exe_path = Sys.iswindows() ? "$outname.exe" : outname
    if !isfile(exe_path)
        bundled_path = joinpath(build_dir, "bin", basename(outname))
        Sys.iswindows() && (bundled_path *= ".exe")
        isfile(bundled_path) && (exe_path = bundled_path)
    end

    n_errors = count(contains("Verifier error"), readlines(logfile))
    return (exe_path, n_errors)
end

function validate_output(exe_path::String, solver_name::String, build_dir::String)
    outfile = joinpath(build_dir, "$(solver_name)_output.csv")
    run(`$(exe_path) $(outfile) 1.5 1.0`)

    lines = readlines(outfile)
    @test length(lines) > 10

    parts = split(lines[1], ",")
    @test length(parts) == 3
    @test parse(Float64, parts[1]) == 0.0
    @test parse(Float64, parts[2]) > 0.0
    @test parse(Float64, parts[3]) > 0.0

    last_parts = split(lines[end], ",")
    @test parse(Float64, last_parts[1]) ≈ 10.0

    return length(lines)
end

# ── Test entry point ───────────────────────────────────────────────────

function run_trim_tests(solvers = SOLVER_ORDER)
    mkpath(LOG_DIR)
    top_build_dir = mktempdir(; cleanup = false)

    for name in solvers
        cfg = SOLVER_CONFIGS[name]
        generate_test_project(name, cfg)
    end

    for name in solvers
        @testset "Trim: $(SOLVER_CONFIGS[name].display_name)" begin
            build_dir = joinpath(top_build_dir, name)
            mkpath(build_dir)

            exe_path, n_errors = compile_solver(name, build_dir)
            @test n_errors == 0
            @test isfile(exe_path)

            if isfile(exe_path) && n_errors == 0
                validate_output(exe_path, name, build_dir)
            end
        end
    end
end

run_trim_tests()

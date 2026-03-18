#!/usr/bin/env julia
#
# Compute which sublibrary test jobs need to run based on changed files.
#
# Reads the dependency graph from lib/*/Project.toml [deps] and [extras]
# sections, identifies internal dependencies by matching dep names against
# known sublibrary directory names, computes the transitive reverse-dependency
# map, then given a list of changed files outputs a GitHub Actions matrix
# include list as JSON.
#
# Each sublibrary can optionally define test groups in test/test_groups.toml:
#
#   [Core]
#   versions = ["lts", "1.11", "1", "pre"]
#
#   [QA]
#   versions = ["1"]
#
#   [GPU]
#   versions = ["1"]
#   runner = ["self-hosted", "Linux", "X64", "gpu"]
#
# Optional fields per group:
#   runner      — string or array of labels (default: "ubuntu-latest")
#   timeout     — integer, job timeout in minutes (default: 120)
#   num_threads — integer, JULIA_NUM_THREADS (default: 1)
#
# If no test/test_groups.toml exists, the default is:
#   Core on ["lts", "1.11", "1", "pre"]
#   QA on ["1"]
#
# Usage:
#   git diff --name-only origin/master...HEAD | julia compute_affected_sublibraries.jl /path/to/repo
#
# Output: JSON array of {group, version, runner, timeout, num_threads} objects
#   for GitHub Actions matrix include.

using TOML

const DEFAULT_TEST_GROUPS = Dict(
    "Core" => ["lts", "1.11", "1", "pre"],
    "QA" => ["1"],
)

function build_dependency_graph(lib_dir::String)
    # Collect all sublibrary names (directories with a Project.toml)
    known_sublibs = Set{String}()
    for entry in readdir(lib_dir)
        if isfile(joinpath(lib_dir, entry, "Project.toml"))
            push!(known_sublibs, entry)
        end
    end

    # Parse each sublibrary's Project.toml for internal deps
    graph = Dict{String, Vector{String}}()
    for pkg in known_sublibs
        toml = TOML.parsefile(joinpath(lib_dir, pkg, "Project.toml"))
        internal_deps = String[]
        for section in ("deps", "extras")
            if haskey(toml, section)
                for dep_name in keys(toml[section])
                    if dep_name in known_sublibs
                        push!(internal_deps, dep_name)
                    end
                end
            end
        end
        graph[pkg] = internal_deps
    end

    return graph
end

function compute_reverse_deps(graph::Dict{String, Vector{String}})
    # Build direct reverse dependency map
    rev = Dict{String, Set{String}}()
    for (pkg, deps) in graph
        for dep in deps
            if !haskey(rev, dep)
                rev[dep] = Set{String}()
            end
            push!(rev[dep], pkg)
        end
    end

    # Compute transitive closure via DFS
    function get_all_rdeps!(visited::Set{String}, pkg::String)
        pkg in visited && return
        push!(visited, pkg)
        for rdep in get(rev, pkg, Set{String}())
            get_all_rdeps!(visited, rdep)
        end
        return
    end

    transitive = Dict{String, Set{String}}()
    for pkg in keys(graph)
        visited = Set{String}()
        get_all_rdeps!(visited, pkg)
        delete!(visited, pkg)  # don't include self
        transitive[pkg] = visited
    end

    return transitive
end

struct TestGroupConfig
    versions::Vector{String}
    runner::Any  # String or Vector{String}
    timeout::Int
    num_threads::Int
end

function load_test_groups(lib_dir::String, pkg::String)
    groups_file = joinpath(lib_dir, pkg, "test", "test_groups.toml")
    if isfile(groups_file)
        toml = TOML.parsefile(groups_file)
        groups = Dict{String, TestGroupConfig}()
        for (name, config) in toml
            versions = convert(Vector{String}, config["versions"])
            runner_raw = get(config, "runner", "ubuntu-latest")
            runner = runner_raw isa Vector ? convert(Vector{String}, runner_raw) : runner_raw::String
            timeout = Int(get(config, "timeout", 120))
            num_threads = Int(get(config, "num_threads", 1))
            groups[name] = TestGroupConfig(versions, runner, timeout, num_threads)
        end
        return groups
    end
    return Dict{String, TestGroupConfig}(
        k => TestGroupConfig(v, "ubuntu-latest", 120, 1) for (k, v) in DEFAULT_TEST_GROUPS
    )
end

function compute_affected(
        changed_files::Vector{String},
        graph::Dict{String, Vector{String}},
        reverse_deps::Dict{String, Set{String}}
    )
    affected = Set{String}()
    for filepath in changed_files
        filepath = strip(filepath)
        isempty(filepath) && continue

        parts = split(filepath, '/')
        if length(parts) >= 2 && parts[1] == "lib" && haskey(graph, String(parts[2]))
            pkg = String(parts[2])
            push!(affected, pkg)
            # Only propagate to reverse deps for src/ or Project.toml changes.
            # Test-only changes don't affect dependents.
            if length(parts) >= 3 &&
                    (parts[3] == "src" || (parts[3] == "Project.toml" && length(parts) == 3))
                union!(affected, get(reverse_deps, pkg, Set{String}()))
            end
        end
    end
    return affected
end

# Entries to exclude from the matrix.
# Each entry is (group, version) where group is the CI GROUP string.
# See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/2977
const EXCLUDES = Set(
    [
        ("OrdinaryDiffEqBDF", "pre"),  # JET resolution fails on pre-release Julia
    ]
)

function build_matrix(affected::Set{String}, lib_dir::String)
    entries = []
    for pkg in sort!(collect(affected))
        groups = load_test_groups(lib_dir, pkg)
        for group_name in sort!(collect(keys(groups)))
            config = groups[group_name]
            ci_group = group_name == "Core" ? pkg : "$(pkg)_$(group_name)"
            for ver in config.versions
                (ci_group, ver) in EXCLUDES && continue
                push!(
                    entries,
                    (;
                        group = ci_group, version = ver, runner = config.runner,
                        timeout = config.timeout, num_threads = config.num_threads,
                    )
                )
            end
        end
    end
    return entries
end

# Minimal JSON serialization (no external dependency needed)
function json_value(v::String)
    return print("\"", v, "\"")
end
function json_value(v::Vector)
    print("[")
    for (j, item) in enumerate(v)
        j > 1 && print(",")
        json_value(item)
    end
    return print("]")
end
function json_value(v::Int)
    return print(v)
end

function print_json(entries)
    print("[")
    for (i, entry) in enumerate(entries)
        i > 1 && print(",")
        print("{\"group\":\"", entry.group, "\",\"version\":\"", entry.version, "\",\"runner\":")
        json_value(entry.runner)
        print(",\"timeout\":", entry.timeout, ",\"num_threads\":", entry.num_threads, "}")
    end
    return println("]")
end

function main()
    if length(ARGS) < 1
        println(stderr, "Usage: julia $(PROGRAM_FILE) <repo_root>")
        exit(1)
    end

    repo_root = ARGS[1]
    lib_dir = joinpath(repo_root, "lib")

    if !isdir(lib_dir)
        println(stderr, "Error: $lib_dir is not a directory")
        exit(1)
    end

    graph = build_dependency_graph(lib_dir)
    reverse_deps = compute_reverse_deps(graph)

    changed_files = split(read(stdin, String), '\n')
    affected = compute_affected(collect(String, changed_files), graph, reverse_deps)

    matrix = build_matrix(affected, lib_dir)
    return print_json(matrix)
end

main()

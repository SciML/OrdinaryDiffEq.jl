#!/usr/bin/env julia
# Run benchmark/benchmarks.jl against a checkout of this monorepo, producing
# a JSON file with `BenchmarkTools.minimum` results.
#
# OrdinaryDiffEq.jl's top-level Project.toml depends on a number of sublibraries
# (OrdinaryDiffEqRosenbrockTableaus, OrdinaryDiffEqExplicitTableaus, …) that
# only live in this repo's `lib/` tree and are not in the General registry.
# They are wired up via `[sources]` in Project.toml. On Julia 1.10 (LTS) Pkg
# does not natively understand `[sources]`, so a plain `Pkg.develop(path=…)`
# of the umbrella package fails to resolve those deps. We therefore parse
# `[sources]` ourselves and `Pkg.develop` each path-source first, then add
# the umbrella package on top.
#
# Usage:
#   julia --startup-file=no run_benchmarks.jl <repo_path> <output_json> [<extra_pkg>…]
#
# Arguments:
#   repo_path     — path to the OrdinaryDiffEq.jl checkout to benchmark.
#   output_json   — file to write JSON results to.
#   extra_pkg…    — additional registered packages to add to the env (e.g.
#                   StableRNGs, StaticArrays, BenchmarkTools is always added).

using Pkg
using TOML

length(ARGS) >= 2 || error("usage: run_benchmarks.jl <repo_path> <output_json> [extra_pkg…]")

const REPO        = abspath(ARGS[1])
const OUTPUT_JSON = abspath(ARGS[2])
const EXTRA_PKGS  = String[ARGS[3:end]...]

isdir(REPO) || error("repo_path is not a directory: $REPO")
isfile(joinpath(REPO, "Project.toml")) ||
    error("no Project.toml at $REPO")
const SCRIPT = joinpath(REPO, "benchmark", "benchmarks.jl")
isfile(SCRIPT) || error("no benchmark script at $SCRIPT")

# Fresh, isolated env so we don't leak state between revs.
Pkg.activate(; temp = true)

# Pre-develop every path-based [sources] entry from the umbrella Project.toml.
# Mirrors AirspeedVelocity's `dev_source_pkgs` so we work on Julia 1.10 too.
project = TOML.parsefile(joinpath(REPO, "Project.toml"))
sources = get(project, "sources", Dict{String, Any}())
path_specs = Pkg.PackageSpec[]
for (name, meta) in sources
    haskey(meta, "path") || continue
    raw = meta["path"]
    p   = isabspath(raw) ? raw : joinpath(REPO, raw)
    push!(path_specs, Pkg.PackageSpec(; name = name, path = p))
end
isempty(path_specs) || Pkg.develop(path_specs)

# The umbrella package itself.
Pkg.develop(Pkg.PackageSpec(; path = REPO))

# Benchmark tooling + script-specific deps.
extras = vcat(["BenchmarkTools"], EXTRA_PKGS)
Pkg.add([Pkg.PackageSpec(; name = p) for p in unique(extras)])

Pkg.instantiate()
Pkg.precompile()

@info "Running benchmark script: $SCRIPT"
include(SCRIPT)  # defines `SUITE::BenchmarkGroup`

@assert @isdefined(SUITE) "benchmark script did not define `SUITE`"

using BenchmarkTools

@info "Tuning benchmarks"
BenchmarkTools.tune!(SUITE; verbose = false)

@info "Running benchmarks"
results = BenchmarkTools.run(SUITE; verbose = true)

# `minimum` collapses each leaf BenchmarkTrial to a single estimate so PR-vs-base
# diffs are stable. `BenchmarkTools.save` writes BenchmarkTools' standard JSON
# wire format, which `BenchmarkTools.load` (and therefore judge/regressions) reads.
mkpath(dirname(OUTPUT_JSON))
BenchmarkTools.save(OUTPUT_JSON, BenchmarkTools.minimum(results))
@info "Wrote results to $OUTPUT_JSON"

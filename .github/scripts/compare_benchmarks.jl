#!/usr/bin/env julia
# Compare two BenchmarkTools JSON result files (as written by
# `BenchmarkTools.save`) and emit a Markdown summary on stdout.
#
# Usage:
#   julia --startup-file=no compare_benchmarks.jl <base.json> <head.json> [<title>]
#
# Output: a Markdown section with one row per benchmark leaf. The "ratio"
# column is `head_time / base_time` — under 1.0 means the PR is faster.

using Pkg
Pkg.activate(; temp = true)
Pkg.add(Pkg.PackageSpec(; name = "BenchmarkTools"))

using BenchmarkTools
using BenchmarkTools: leaves
using Printf

length(ARGS) >= 2 || error("usage: compare_benchmarks.jl <base.json> <head.json> [<title>]")

const BASE_PATH = ARGS[1]
const HEAD_PATH = ARGS[2]
const TITLE     = length(ARGS) >= 3 ? ARGS[3] : "Benchmark results"

isfile(BASE_PATH) || error("base results not found: $BASE_PATH")
isfile(HEAD_PATH) || error("head results not found: $HEAD_PATH")

base = BenchmarkTools.load(BASE_PATH)[1]::BenchmarkGroup
head = BenchmarkTools.load(HEAD_PATH)[1]::BenchmarkGroup

# Format a duration in human-friendly units.
function fmt_time(ns::Real)
    ns < 1e3   && return @sprintf("%.2f ns", ns)
    ns < 1e6   && return @sprintf("%.2f μs", ns / 1e3)
    ns < 1e9   && return @sprintf("%.2f ms", ns / 1e6)
    return @sprintf("%.2f s", ns / 1e9)
end
fmt_mem(b::Real) = b < 1024 ? @sprintf("%d B", b) :
    b < 1024^2 ? @sprintf("%.1f KiB", b / 1024) :
    b < 1024^3 ? @sprintf("%.2f MiB", b / 1024^2) :
                  @sprintf("%.2f GiB", b / 1024^3)

# `leaves` returns Vector{Tuple{Vector{Any}, Trial-or-Estimate}}. Build a
# {key_path => estimate} dict for each side so we can join.
to_map(g) = Dict(join(string.(k), " / ") => v for (k, v) in leaves(g))
base_map = to_map(base)
head_map = to_map(head)
all_keys = sort!(collect(union(keys(base_map), keys(head_map))))

println("## ", TITLE)
println()
println("| Benchmark | base time | PR time | ratio | base mem | PR mem |")
println("|---|---:|---:|---:|---:|---:|")

regressed = String[]
for k in all_keys
    b = get(base_map, k, nothing)
    h = get(head_map, k, nothing)
    bt = b === nothing ? NaN : time(b)
    ht = h === nothing ? NaN : time(h)
    bm = b === nothing ? 0   : memory(b)
    hm = h === nothing ? 0   : memory(h)
    ratio = (isnan(bt) || isnan(ht) || bt == 0) ? NaN : ht / bt
    ratio_cell = isnan(ratio) ? "—" : @sprintf("%.2fx", ratio)
    if !isnan(ratio) && ratio > 1.10
        ratio_cell *= " ⚠️"
        push!(regressed, k)
    elseif !isnan(ratio) && ratio < 0.90
        ratio_cell *= " 🚀"
    end
    println("| `", k, "` | ",
        isnan(bt) ? "—" : fmt_time(bt), " | ",
        isnan(ht) ? "—" : fmt_time(ht), " | ",
        ratio_cell, " | ",
        b === nothing ? "—" : fmt_mem(bm), " | ",
        h === nothing ? "—" : fmt_mem(hm), " |")
end

println()
if isempty(regressed)
    println("_No benchmarks regressed by more than 10%._")
else
    println("**", length(regressed), " benchmark(s) regressed by more than 10%:**")
    for k in regressed
        println("- `", k, "`")
    end
end

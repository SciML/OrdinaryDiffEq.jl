using Documenter
using StochasticDiffEq
# StochasticDiffEq @reexports its solver subpackages, but the docstrings the @docs
# blocks reference live in those subpackages, so they must be in `modules` for
# Documenter to find them.
using StochasticDiffEqCore, StochasticDiffEqHighOrder, StochasticDiffEqIIF,
    StochasticDiffEqImplicit, StochasticDiffEqLeaping, StochasticDiffEqLowOrder,
    StochasticDiffEqMilstein, StochasticDiffEqROCK, StochasticDiffEqRODE,
    StochasticDiffEqWeak

cp(
    joinpath(@__DIR__, "Project.toml"), joinpath(@__DIR__, "src", "assets", "Project.toml"),
    force = true
)

# Keep pages.jl separate so DiffEqDocs.jl can include it when aggregating these docs.
include("pages.jl")

makedocs(;
    sitename = "StochasticDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    checkdocs = :exports,
    modules = [
        StochasticDiffEq, StochasticDiffEqCore, StochasticDiffEqHighOrder,
        StochasticDiffEqIIF, StochasticDiffEqImplicit, StochasticDiffEqLeaping,
        StochasticDiffEqLowOrder, StochasticDiffEqMilstein, StochasticDiffEqROCK,
        StochasticDiffEqRODE, StochasticDiffEqWeak,
    ],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/StochasticDiffEq/stable/"
    ),
    pages,
)

# Note: these pages are aggregated into the unified SciML docs by DiffEqDocs.jl
# (which copies docs/src + docs/pages.jl); deployment happens there.

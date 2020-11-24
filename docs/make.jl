using Documenter, OrdinaryDiffEq

makedocs(
    sitename="OrdinaryDiffEq.jl",
    authors="Chris Rackauckas et al.",
    clean=true,
    doctest=false,
    modules=[OrdinaryDiffEq],

    format=Documenter.HTML(assets=["assets/favicon.ico"],
                           canonical="https://ordinarydiffeq.sciml.ai/stable/"),

    pages=[
        "OrdinaryDiffEq.jl: ODE solvers and utilities" => "index.md",
        "Usage" => "usage.md"
    ]
)

deploydocs(
    repo="github.com/SciML/OrdinaryDiffEq.jl";
    push_preview=true
)

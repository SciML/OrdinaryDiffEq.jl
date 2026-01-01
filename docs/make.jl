using Documenter, OrdinaryDiffEq

# Handle both running from repo root and from docs directory
manifest_src = isfile("./docs/Manifest.toml") ? "./docs/Manifest.toml" : "Manifest.toml"
project_src = isfile("./docs/Project.toml") ? "./docs/Project.toml" : "Project.toml"
manifest_dst = isfile("./docs/src/assets/Manifest.toml") || isdir("./docs/src/assets") ? "./docs/src/assets/Manifest.toml" : "src/assets/Manifest.toml"
project_dst = isfile("./docs/src/assets/Project.toml") || isdir("./docs/src/assets") ? "./docs/src/assets/Project.toml" : "src/assets/Project.toml"

cp(manifest_src, manifest_dst, force = true)
cp(project_src, project_dst, force = true)

# Keep pages.jl separate for the DiffEqDocs.jl build
include("pages.jl")

makedocs(sitename = "OrdinaryDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    modules = [OrdinaryDiffEq,
        OrdinaryDiffEq.OrdinaryDiffEqAdamsBashforthMoulton,
        OrdinaryDiffEq.OrdinaryDiffEqBDF,
        OrdinaryDiffEq.OrdinaryDiffEqDefault,
        OrdinaryDiffEq.OrdinaryDiffEqExplicitRK,
        OrdinaryDiffEq.OrdinaryDiffEqExponentialRK,
        OrdinaryDiffEq.OrdinaryDiffEqExtrapolation,
        OrdinaryDiffEq.OrdinaryDiffEqFeagin,
        OrdinaryDiffEq.OrdinaryDiffEqFIRK,
        OrdinaryDiffEq.OrdinaryDiffEqHighOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqIMEXMultistep,
        OrdinaryDiffEq.OrdinaryDiffEqLinear,
        OrdinaryDiffEq.OrdinaryDiffEqLowOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqLowStorageRK,
        OrdinaryDiffEq.OrdinaryDiffEqNordsieck,
        OrdinaryDiffEq.OrdinaryDiffEqPDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqPRK,
        OrdinaryDiffEq.OrdinaryDiffEqQPRK,
        OrdinaryDiffEq.OrdinaryDiffEqRKN,
        OrdinaryDiffEq.OrdinaryDiffEqRosenbrock,
        OrdinaryDiffEq.OrdinaryDiffEqSDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqSSPRK,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEq.OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEq.OrdinaryDiffEqTsit5,
        OrdinaryDiffEq.OrdinaryDiffEqVerner
    ],
    warnonly = [:docs_block, :missing_docs, :eval_block],
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://ordinarydiffeq.sciml.ai/stable/",
        size_threshold_ignore = [joinpath("semiimplicit", "Rosenbrock.md"),
            joinpath("massmatrixdae", "Rosenbrock.md")]),
    pages = pages)

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true)

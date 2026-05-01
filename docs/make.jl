using Documenter, OrdinaryDiffEq, DiffEqDevTools
using OrdinaryDiffEqAMF
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqBDF
using OrdinaryDiffEqDefault
using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExponentialRK
using OrdinaryDiffEqExtrapolation
using OrdinaryDiffEqFeagin
using OrdinaryDiffEqFIRK
using OrdinaryDiffEqHighOrderRK
using OrdinaryDiffEqIMEXMultistep
using OrdinaryDiffEqLinear
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqNordsieck
using OrdinaryDiffEqPDIRK
using OrdinaryDiffEqPRK
using OrdinaryDiffEqQPRK
using OrdinaryDiffEqRKN
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqStabilizedIRK
using OrdinaryDiffEqStabilizedRK
using OrdinaryDiffEqSymplecticRK
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner

cp(joinpath(@__DIR__, "Manifest.toml"), joinpath(@__DIR__, "src", "assets", "Manifest.toml"), force = true)
cp(joinpath(@__DIR__, "Project.toml"), joinpath(@__DIR__, "src", "assets", "Project.toml"), force = true)

# Keep pages.jl separate for the DiffEqDocs.jl build
include("pages.jl")

makedocs(
    sitename = "OrdinaryDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    modules = [
        OrdinaryDiffEq,
        OrdinaryDiffEqAdamsBashforthMoulton,
        OrdinaryDiffEqBDF,
        OrdinaryDiffEqDefault,
        OrdinaryDiffEqExplicitRK,
        OrdinaryDiffEqExponentialRK,
        OrdinaryDiffEqExtrapolation,
        OrdinaryDiffEqFeagin,
        OrdinaryDiffEqFIRK,
        OrdinaryDiffEqHighOrderRK,
        OrdinaryDiffEqIMEXMultistep,
        OrdinaryDiffEqLinear,
        OrdinaryDiffEqLowOrderRK,
        OrdinaryDiffEqLowStorageRK,
        OrdinaryDiffEqNordsieck,
        OrdinaryDiffEqPDIRK,
        OrdinaryDiffEqPRK,
        OrdinaryDiffEqQPRK,
        OrdinaryDiffEqRKN,
        OrdinaryDiffEqRosenbrock,
        OrdinaryDiffEqSDIRK,
        OrdinaryDiffEqSSPRK,
        OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEqTsit5,
        OrdinaryDiffEqVerner,
        OrdinaryDiffEqAMF,
        DiffEqDevTools,
    ],
    linkcheck_ignore = [r"https://github.com/JuliaDiff/ForwardDiff.jl"],
    warnonly = [:docs_block, :missing_docs, :eval_block],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://ordinarydiffeq.sciml.ai/stable/",
        size_threshold_ignore = [
            joinpath("semiimplicit", "Rosenbrock.md"),
            joinpath("massmatrixdae", "Rosenbrock.md"),
        ]
    ),
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true
)

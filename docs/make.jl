using Documenter, OrdinaryDiffEq

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

makedocs(sitename = "OrdinaryDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    modules = [OrdinaryDiffEq,
        OrdinaryDiffEq.OrdinaryDiffEqExtrapolation,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEq.OrdinaryDiffEqLowStorageRK,
        OrdinaryDiffEq.OrdinaryDiffEqSSPRK,
        OrdinaryDiffEq.OrdinaryDiffEqFeagin,
        OrdinaryDiffEq.OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEq.OrdinaryDiffEqRKN,
        OrdinaryDiffEq.OrdinaryDiffEqTsit5,
        OrdinaryDiffEq.OrdinaryDiffEqLowOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqHighOrderRK,
        OrdinaryDiffEq.OrdinaryDiffEqVerner,
        OrdinaryDiffEq.OrdinaryDiffEqAdamsBashforthMoulton,
        OrdinaryDiffEq.OrdinaryDiffEqNordsieck,
        OrdinaryDiffEq.OrdinaryDiffEqPRK,
        OrdinaryDiffEq.OrdinaryDiffEqRKN,
        OrdinaryDiffEq.OrdinaryDiffEqSDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqBDF,
        OrdinaryDiffEq.OrdinaryDiffEqDefault,
        OrdinaryDiffEq.OrdinaryDiffEqFIRK,
        OrdinaryDiffEq.OrdinaryDiffEqPDIRK,
        OrdinaryDiffEq.OrdinaryDiffEqRosenbrock,
        OrdinaryDiffEq.OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEq.OrdinaryDiffEqBDF],
    warnonly = [:docs_block, :missing_docs, :eval_block],
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://ordinarydiffeq.sciml.ai/stable/"),
    pages = [
        "OrdinaryDiffEq.jl: ODE solvers and utilities" => "index.md",
        "Usage" => "usage.md",
        "Explicit Solvers" => [
            "explicit/Tsit5.md",
            "explicit/LowOrderRK.md",
            "explicit/HighOrderRK.md",
            "explicit/Verner.md",
            "explicit/Feagin.md",
            "explicit/LowStorageRK.md",
            "explicit/SSPRK.md",
            "explicit/AdamsBashforthMoulton.md",
            "explicit/Nordsieck.md",
            "explicit/RKN.md",
            "explicit/SymplecticRK.md",
            "explicit/PRK.md",
            "explicit/Extrapolation.md",
        ],
        "Implicit Solvers" => [
            "implicit/SDIRK.md",
            "implicit/FIRK.md",
            "implicit/PDIRK.md",
            "implicit/Rosenbrock.md",
            "implicit/StabalizedRK.md",
            "implicit/BDF.md",
        ],
        "IMEX Solvers" => [
            "imex/imex_multistep.md",
            "imex/imex_sdirk.md"
        ],
        "Semilinear ODE Solvers" => [
            "semilinear/exponential_rk.md",
            "semilinear/magnus.md"
        ],
        "Misc Solvers" => [
            "misc.md"
        ]
    ])

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true)

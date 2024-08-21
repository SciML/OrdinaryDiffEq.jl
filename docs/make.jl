using Documenter, OrdinaryDiffEq

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

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
            "explicit/QPRK.md",
            "explicit/Extrapolation.md"
        ],
        "Implicit Solvers" => [
            "implicit/SDIRK.md",
            "implicit/FIRK.md",
            "implicit/PDIRK.md",
            "implicit/Rosenbrock.md",
            "implicit/StabalizedRK.md",
            "implicit/StabalizedIRK.md",
            "implicit/BDF.md",
            "implicit/Extrapolation.md"
        ],
        "IMEX Solvers" => [
            "imex/IMEXMultistep.md"
        ],
        "Semilinear ODE Solvers" => [
            "semilinear/ExponentialRK.md",
            "semilinear/Linear.md"
        ],
        "Misc Solvers" => [
            "misc.md"
        ]
    ])

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true)

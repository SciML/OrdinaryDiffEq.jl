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
        canonical = "https://ordinarydiffeq.sciml.ai/stable/",
        size_threshold_ignore = [joinpath("implicit", "Rosenbrock.md"),
                            joinpath("massmatrixdae", "Rosenbrock.md")],
    pages = [
        "OrdinaryDiffEq.jl: ODE solvers and utilities" => "index.md",
        "Usage" => "usage.md",
        "Explicit Solvers" => [
            "explicit/Tsit5.md",
            "explicit/Verner.md",
            "explicit/AdamsBashforthMoulton.md",
            "explicit/LowStorageRK.md",
            "explicit/SSPRK.md",
            "explicit/LowOrderRK.md",
            "explicit/HighOrderRK.md",
            "explicit/Feagin.md",
            "explicit/PRK.md",
            "explicit/QPRK.md",
            "explicit/Extrapolation.md"
        ],
        "Semi-Implicit Solvers" => [
            "semiimplicit/Rosenbrock.md",
            "semiimplicit/StabalizedRK.md",
            "semiimplicit/ExponentialRK.md",
        ],
        "Implicit Solvers" => [
            "implicit/SDIRK.md",
            "implicit/FIRK.md",
            "implicit/BDF.md",
            "implicit/Extrapolation.md",
            "implicit/PDIRK.md",
            "implicit/Nordsieck.md",
            ],
        "IMEX Solvers" => [
            "imex/IMEXMultistep.md",
            "imex/StabalizedIRK.md",
            "imex/IMEXBDF.md"
        ],
        "Dynamical ODE Explicit Solvers" => [
            "dynamicalodeexplicit/RKN.md",
            "dynamicalodeexplicit/SymplecticRK.md",
        ],
        "Semilinear ODE Solvers" => [
            "semilinear/ExponentialRK.md",
            "semilinear/Linear.md"
        ],
        "Mass Matrix DAE Solvers" => [
            "massmatrixdae/Rosenbrock.md",
            "massmatrixdae/BDF.md",
        ],
        "Fully Implicit DAE Solvers" => [
            "fullyimplicitdae/BDF.md",
        ],
        "Misc Solvers" => [
            "misc.md"
        ]
    ])

deploydocs(repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true)

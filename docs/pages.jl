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
        "explicit/TaylorSeries.md",
        "explicit/Extrapolation.md",
    ],
    "Semi-Implicit Solvers" => [
        "semiimplicit/Rosenbrock.md",
        "semiimplicit/AMF.md",
        "semiimplicit/StabilizedRK.md",
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
        "imex/StabilizedIRK.md",
        "imex/IMEXBDF.md",
    ],
    "Dynamical ODE Explicit Solvers" => [
        "dynamicalodeexplicit/RKN.md",
        "dynamicalodeexplicit/SymplecticRK.md",
    ],
    "Semilinear ODE Solvers" => [
        "semilinear/ExponentialRK.md",
        "semilinear/Linear.md",
    ],
    "Mass Matrix DAE Solvers" => [
        "massmatrixdae/Rosenbrock.md",
        "massmatrixdae/BDF.md",
    ],
    "Fully Implicit DAE Solvers" => [
        "fullyimplicitdae/BDF.md",
    ],
    "Misc Solvers" => [
        "misc.md",
    ],
    "Developer Documentation" => [
        "devtools/index.md",
        "Contributor Guide" => [
            "devtools/contributing/ecosystem_overview.md",
            "devtools/contributing/adding_packages.md",
            "devtools/contributing/adding_algorithms.md",
            "devtools/contributing/defining_problems.md",
            "devtools/contributing/diffeq_internals.md",
            "devtools/contributing/type_traits.md",
        ],
        "Algorithm Development Tools" => [
            "devtools/alg_dev/test_problems.md",
            "devtools/alg_dev/convergence.md",
            "devtools/alg_dev/benchmarks.md",
        ],
        "Internal Documentation" => [
            "devtools/internals/notes_on_algorithms.md",
            "devtools/internals/tableaus.md",
        ],
    ],
]

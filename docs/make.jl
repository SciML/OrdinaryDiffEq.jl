using Documenter, OrdinaryDiffEq, DiffEqDevTools
import ADTypes
import CommonSolve
using DelayDiffEq
import SciMLBase, SciMLLogging, SciMLOperators
using DiffEqBase
using OrdinaryDiffEqCore
# Bring controller API symbols into Main so unqualified @ref links in
# docs/src/api/controllers.md resolve. These are not exported by
# OrdinaryDiffEqCore but are documented public API.
using OrdinaryDiffEqCore: default_controller, resolve_basic,
    get_EEst, set_EEst!, CompositeController
using OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqFunctionMap
using ImplicitDiscreteSolve
using StochasticDiffEqLevyArea
using StochasticDiffEqWeak
using StochasticDiffEqCore
using StochasticDiffEqHighOrder
using StochasticDiffEqIIF
using StochasticDiffEqImplicit
using StochasticDiffEqLeaping
using StochasticDiffEqLowOrder
using StochasticDiffEqMilstein
using StochasticDiffEqROCK
using StochasticDiffEqRODE

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
using OrdinaryDiffEqNewmark
using OrdinaryDiffEqRKN
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqRosenbrockTableaus
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqStabilizedIRK
using OrdinaryDiffEqStabilizedRK
using OrdinaryDiffEqSymplecticRK
import OrdinaryDiffEqTaylorSeries
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner
using GlobalDiffEq

cp(joinpath(@__DIR__, "Manifest.toml"), joinpath(@__DIR__, "src", "assets", "Manifest.toml"), force = true)
cp(joinpath(@__DIR__, "Project.toml"), joinpath(@__DIR__, "src", "assets", "Project.toml"), force = true)

# Keep pages.jl separate for the DiffEqDocs.jl build
include("pages.jl")

makedocs(;
    sitename = "OrdinaryDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    modules = [
        OrdinaryDiffEq,
        DiffEqBase,
        OrdinaryDiffEqCore,
        OrdinaryDiffEqDifferentiation,
        OrdinaryDiffEqNonlinearSolve,
        OrdinaryDiffEqFunctionMap,
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
        OrdinaryDiffEqNewmark,
        OrdinaryDiffEqRKN,
        OrdinaryDiffEqRosenbrock,
        OrdinaryDiffEqRosenbrockTableaus,
        OrdinaryDiffEqSDIRK,
        OrdinaryDiffEqSSPRK,
        OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEqTaylorSeries,
        OrdinaryDiffEqTsit5,
        OrdinaryDiffEqVerner,
        OrdinaryDiffEqAMF,
        ImplicitDiscreteSolve,
        StochasticDiffEqLevyArea,
        StochasticDiffEqCore,
        StochasticDiffEqHighOrder,
        StochasticDiffEqIIF,
        StochasticDiffEqImplicit,
        StochasticDiffEqLeaping,
        StochasticDiffEqLowOrder,
        StochasticDiffEqMilstein,
        StochasticDiffEqROCK,
        StochasticDiffEqRODE,
        StochasticDiffEqWeak,
        DiffEqDevTools,
        DelayDiffEq,
        GlobalDiffEq,
    ],
    checkdocs = :public,
    linkcheck_ignore = [r"https://github.com/JuliaDiff/ForwardDiff.jl"],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://ordinarydiffeq.sciml.ai/stable/",
        size_threshold_ignore = [
            joinpath("semiimplicit", "Rosenbrock.md"),
            joinpath("massmatrixdae", "Rosenbrock.md"),
            joinpath("devtools", "internals", "public_api.md"),
        ]
    ),
    pages
)

deploydocs(
    repo = "github.com/SciML/OrdinaryDiffEq.jl";
    push_preview = true
)

using Documenter, OrdinaryDiffEq, DiffEqDevTools
import SciMLBase, SciMLLogging
using DiffEqBase
using OrdinaryDiffEqCore
# Bring controller API symbols into Main so unqualified @ref links in
# docs/src/api/controllers.md resolve. These are not exported by
# OrdinaryDiffEqCore but are documented public API.
using OrdinaryDiffEqCore: default_controller, resolve_basic,
    get_EEst, set_EEst!, CompositeController
using OrdinaryDiffEqNonlinearSolve
using ImplicitDiscreteSolve

# Register the re-exported bindings with the umbrella module for Documenter.
for (name, owner) in (
        (:solve, "CommonSolve.solve"),
        (:solve!, "CommonSolve.solve!"),
        (:init, "CommonSolve.init"),
        (:step!, "CommonSolve.step!"),
        (:ODEProblem, "SciMLBase.ODEProblem"),
        (:ODEFunction, "SciMLBase.ODEFunction"),
        (:ODESolution, "SciMLBase.ODESolution"),
        (:SplitODEProblem, "SciMLBase.SplitODEProblem"),
        (:SplitFunction, "SciMLBase.SplitFunction"),
        (:SecondOrderODEProblem, "SciMLBase.SecondOrderODEProblem"),
        (:DynamicalODEProblem, "SciMLBase.DynamicalODEProblem"),
        (:DAEProblem, "SciMLBase.DAEProblem"),
        (:DAEFunction, "SciMLBase.DAEFunction"),
        (:DAESolution, "SciMLBase.DAESolution"),
        (:EnsembleProblem, "SciMLBase.EnsembleProblem"),
        (:CallbackSet, "SciMLBase.CallbackSet"),
        (:ContinuousCallback, "SciMLBase.ContinuousCallback"),
        (:DiscreteCallback, "SciMLBase.DiscreteCallback"),
        (:VectorContinuousCallback, "SciMLBase.VectorContinuousCallback"),
        (:ODEAliasSpecifier, "SciMLBase.ODEAliasSpecifier"),
        (:add_tstop!, "SciMLBase.add_tstop!"),
        (:derivative_discontinuity!, "SciMLBase.derivative_discontinuity!"),
        (:reinit!, "SciMLBase.reinit!"),
        (:remake, "SciMLBase.remake"),
        (:set_proposed_dt!, "SciMLBase.set_proposed_dt!"),
        (:successful_retcode, "SciMLBase.successful_retcode"),
        (:AutoFiniteDiff, "ADTypes.AutoFiniteDiff"),
        (:AutoForwardDiff, "ADTypes.AutoForwardDiff"),
        (:AutoSparse, "ADTypes.AutoSparse"),
        (:DefaultODEAlgorithm, "OrdinaryDiffEqDefault.DefaultODEAlgorithm"),
    )
    doc = "Re-export of `$owner`."
    Core.eval(OrdinaryDiffEq, :(@doc $doc $name))
end

@eval SciMLBase @doc "The common interfaces and types used throughout the SciML ecosystem." SciMLBase
@eval SciMLLogging @doc "The common logging interface used by SciML solvers." SciMLLogging
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
        SciMLBase,
        SciMLBase.ReturnCode,
        SciMLLogging,
        DiffEqBase,
        OrdinaryDiffEqCore,
        OrdinaryDiffEqNonlinearSolve,
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
        OrdinaryDiffEqSDIRK,
        OrdinaryDiffEqSSPRK,
        OrdinaryDiffEqStabilizedIRK,
        OrdinaryDiffEqStabilizedRK,
        OrdinaryDiffEqSymplecticRK,
        OrdinaryDiffEqTsit5,
        OrdinaryDiffEqVerner,
        OrdinaryDiffEqAMF,
        ImplicitDiscreteSolve,
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

module DelayDiffEq

import Reexport: @reexport, Reexport
@reexport using SciMLBase
import OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqDifferentiation,
    OrdinaryDiffEqRosenbrock
import OrdinaryDiffEqDefault: OrdinaryDiffEqDefault
import OrdinaryDiffEqFunctionMap: OrdinaryDiffEqFunctionMap

using DataStructures: BinaryMinHeap
using LinearAlgebra: opnorm, I
using Logging: @logmsg
using RecursiveArrayTools: copyat_or_push!, recursivecopy, recursivecopy!,
    recursive_bottom_eltype, recursive_unitless_bottom_eltype,
    recursive_unitless_eltype
using ForwardDiff: ForwardDiff

import ArrayInterface
import SimpleNonlinearSolve
import SymbolicIndexingInterface as SII

using SciMLBase: AbstractDDEAlgorithm, AbstractDDEIntegrator, AbstractODEIntegrator,
    DEIntegrator
using SciMLBase: AbstractSDDEProblem, SDDEProblem, AbstractSDDEAlgorithm

using Base: deleteat!
import FastBroadcast: @..

using OrdinaryDiffEqNonlinearSolve: NLAnderson, NLFunctional
using OrdinaryDiffEqCore: AbstractNLSolverCache, SlowConvergence,
    alg_extrapolates, alg_maximum_order, initialize!, DEVerbosity
using OrdinaryDiffEqCore: StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqRODEAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache
using OrdinaryDiffEqRosenbrock: RosenbrockMutableCache
using OrdinaryDiffEqFunctionMap: FunctionMap
# using OrdinaryDiffEqDifferentiation: resize_grad_config!, resize_jac_config!

using DiffEqBase: is_diagonal_noise

# Explicit imports for functions
using OrdinaryDiffEqCore: AutoSwitch, CompositeAlgorithm
using OrdinaryDiffEqDefault: DefaultODEAlgorithm
using SciMLBase: CallbackSet, DAEProblem, DDEProblem, DESolution, ODEProblem, ReturnCode,
    VectorContinuousCallback, addat!, addat_non_user_cache!,
    deleteat_non_user_cache!,
    du_cache, full_cache, get_tmp_cache, isinplace,
    reeval_internals_due_to_modification!,
    reinit!, resize_non_user_cache!, savevalues!, u_cache, user_cache,
    step!, terminate!, u_modified!, get_proposed_dt, set_proposed_dt!,
    auto_dt_reset!,
    add_tstop!, add_saveat!, get_du, get_du!, addsteps!,
    change_t_via_interpolation!, isadaptive
using DiffEqBase: initialize!
import DiffEqBase
using SciMLLogging: AbstractVerbosityPreset, None, @SciMLMessage

import SciMLBase

const SDEAlgUnion = Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}

# Internal hook functions for SDDE support. These are overloaded by the
# StochasticDiffEqCore extension to provide actual implementations.
# Calling them without the extension loaded gives a clear error.
function _sde_alg_cache end
function _create_sdde_noise end

export Discontinuity, MethodOfSteps

include("discontinuity_type.jl")
include("functionwrapper.jl")

include("integrators/type.jl")
include("integrators/utils.jl")
include("integrators/interface.jl")

include("fpsolve/type.jl")
include("fpsolve/fpsolve.jl")
include("fpsolve/utils.jl")
include("fpsolve/functional.jl")
include("cache_utils.jl")
include("interpolants.jl")
include("history_function.jl")
include("algorithms.jl")
include("track.jl")
include("alg_utils.jl")
include("solve.jl")
include("utils.jl")

# Default solver for DDEProblems (no algorithm specified)
function SciMLBase.__solve(prob::DDEProblem; kwargs...)
    return SciMLBase.solve(prob, MethodOfSteps(DefaultODEAlgorithm()); kwargs...)
end
function SciMLBase.__init(prob::DDEProblem; kwargs...)
    return DiffEqBase.init(prob, MethodOfSteps(DefaultODEAlgorithm()); kwargs...)
end

# Auto-wrap bare ODE/SDE algorithms in MethodOfSteps for DDE/SDDE problems.
# Allows: solve(DDEProblem(...), Tsit5()) instead of MethodOfSteps(Tsit5())
# and:    solve(SDDEProblem(...), EM()) instead of MethodOfSteps(EM())
const OrdinaryDiffEqAlgorithm = OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm
function SciMLBase.__solve(
        prob::SciMLBase.AbstractDDEProblem,
        alg::OrdinaryDiffEqAlgorithm, args...; kwargs...
    )
    return SciMLBase.__solve(prob, MethodOfSteps(alg), args...; kwargs...)
end
function SciMLBase.__init(
        prob::SciMLBase.AbstractDDEProblem,
        alg::OrdinaryDiffEqAlgorithm, args...; kwargs...
    )
    return SciMLBase.__init(prob, MethodOfSteps(alg), args...; kwargs...)
end
function DiffEqBase.check_prob_alg_pairing(prob::DDEProblem, alg::OrdinaryDiffEqAlgorithm)
    return nothing
end

function SciMLBase.__solve(
        prob::AbstractSDDEProblem,
        alg::SDEAlgUnion, args...; kwargs...
    )
    return SciMLBase.__solve(prob, MethodOfSteps(alg), args...; kwargs...)
end
function SciMLBase.__init(
        prob::AbstractSDDEProblem,
        alg::SDEAlgUnion, args...; kwargs...
    )
    return SciMLBase.__init(prob, MethodOfSteps(alg), args...; kwargs...)
end
function DiffEqBase.check_prob_alg_pairing(prob::SDDEProblem, alg::SDEAlgUnion)
    return nothing
end

end # module

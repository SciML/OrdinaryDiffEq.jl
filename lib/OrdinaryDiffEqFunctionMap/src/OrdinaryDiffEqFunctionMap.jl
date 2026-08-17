module OrdinaryDiffEqFunctionMap

import OrdinaryDiffEqCore: isfsal, beta2_default, beta1_default, OrdinaryDiffEqAlgorithm,
    perform_step!, unwrap_alg,
    OrdinaryDiffEqMutableCache,
    alg_cache, @cache, _ode_addsteps!, _ode_interpolant!,
    _ode_interpolant, get_fsalfirstlast,
    OrdinaryDiffEqConstantCache, dt_required,
    isdiscretecache, isdiscretealg,
    trivial_limiter!
using DiffEqBase: DiffEqBase
import DiffEqBase: initialize!
import RecursiveArrayTools: recursivecopy!
import FastBroadcast: @..
import MuladdMacro: @muladd
import OrdinaryDiffEqCore

import SciMLBase
import SciMLBase: alg_order

using Reexport: Reexport, @reexport
@reexport using SciMLBase: ODEProblem, ODEFunction, solve, init, solve!, step!, remake,
    reinit!, ReturnCode, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    CallbackSet, terminate!, add_tstop!, derivative_discontinuity!, set_proposed_dt!,
    successful_retcode, ODEAliasSpecifier, DiscreteProblem, DiscreteFunction

include("algorithms.jl")
include("alg_utils.jl")
include("functionmap_caches.jl")
include("interp_func.jl")
include("interpolants.jl")
include("functionmap_perform_step.jl")
include("fixed_timestep_perform_step.jl")

# Default algorithm for DiscreteProblem
function SciMLBase.__solve(
        prob::SciMLBase.DiscreteProblem, ::Nothing, args...;
        kwargs...
    )
    return SciMLBase.__solve(prob, FunctionMap(), args...; kwargs...)
end

function SciMLBase.__init(
        prob::SciMLBase.DiscreteProblem, ::Nothing, args...;
        kwargs...
    )
    return SciMLBase.__init(prob, FunctionMap(), args...; kwargs...)
end

export FunctionMap

end

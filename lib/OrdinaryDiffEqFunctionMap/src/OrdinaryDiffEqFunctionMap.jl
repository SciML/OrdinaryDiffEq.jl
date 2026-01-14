module OrdinaryDiffEqFunctionMap

import OrdinaryDiffEqCore: ODEVerbosity, isfsal, beta2_default, beta1_default, OrdinaryDiffEqAlgorithm,
    initialize!, perform_step!, unwrap_alg,
    OrdinaryDiffEqMutableCache,
    alg_cache, @cache, _ode_addsteps!, _ode_interpolant!,
    _ode_interpolant, get_fsalfirstlast,
    alg_order, OrdinaryDiffEqConstantCache, dt_required,
    isdiscretecache, isdiscretealg, full_cache
using DiffEqBase
import RecursiveArrayTools: recursivecopy!
import FastBroadcast: @..
import MuladdMacro: @muladd
import Static: False
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

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

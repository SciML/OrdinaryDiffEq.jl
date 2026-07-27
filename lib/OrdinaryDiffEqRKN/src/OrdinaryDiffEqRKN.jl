module OrdinaryDiffEqRKN

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    perform_step!,
    alg_extrapolates,
    OrdinaryDiffEqAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptivePartitionedAlgorithm,
    OrdinaryDiffEqPartitionedAlgorithm,
    CompiledFloats,
    alg_cache, @cache, full_cache,
    constvalue, _ode_interpolant,
    get_fsalfirstlast,
    _ode_interpolant!,
    generic_solver_docstring
import SciMLBase: alg_order, @def
using SciMLBase: SciMLBase
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!,
    @tight_loop_macros
using FastBroadcast: FastBroadcast, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: RecursiveArrayTools, ArrayPartition, recursivefill!
import OrdinaryDiffEqCore

using Reexport: Reexport, @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("rkn_tableaus.jl")
include("rkn_caches.jl")
include("interp_func.jl")
include("interpolants.jl")
include("rkn_perform_step.jl")

export Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
    Nystrom5VelocityIndependent,
    IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
    RKN4

end

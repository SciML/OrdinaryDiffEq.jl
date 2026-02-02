module ImplicitDiscreteSolve

using SciMLBase
using NonlinearSolveBase
using NonlinearSolveFirstOrder
using SymbolicIndexingInterface: parameter_symbols
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, alg_cache, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache, get_fsalfirstlast, isfsal,
    initialize!, perform_step!, isdiscretecache, isdiscretealg,
    alg_order, beta2_default, beta1_default, dt_required,
    _initialize_dae!, DefaultInit, BrownFullBasicInit, OverrideInit,
    @muladd, @.., _unwrap_val, OrdinaryDiffEqCore, isadaptive,
    AbstractController

import OrdinaryDiffEqCore: AbstractControllerCache

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("cache.jl")
include("solve.jl")
include("alg_utils.jl")
include("controller.jl")

export IDSolve

end

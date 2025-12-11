module ImplicitDiscreteSolve

using SciMLBase
using SimpleNonlinearSolve
using SymbolicIndexingInterface: parameter_symbols
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, alg_cache, OrdinaryDiffEqMutableCache,
                           OrdinaryDiffEqConstantCache, get_fsalfirstlast, isfsal,
                           initialize!, perform_step!, isdiscretecache, isdiscretealg,
                           alg_order, beta2_default, beta1_default, dt_required,
                           _initialize_dae!, DefaultInit, BrownFullBasicInit, OverrideInit

using Reexport
@reexport using SciMLBase

"""
    IDSolve()

Simple solver for `ImplicitDiscreteSystems`. Uses `SimpleNewtonRaphson` to solve for the next state at every timestep.
"""
struct IDSolve <: OrdinaryDiffEqAlgorithm end

include("cache.jl")
include("solve.jl")
include("alg_utils.jl")

export IDSolve

end

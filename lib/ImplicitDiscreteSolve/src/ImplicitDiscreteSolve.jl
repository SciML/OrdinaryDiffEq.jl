module ImplicitDiscreteSolve

using SciMLBase: SciMLBase, ImplicitDiscreteProblem, NonlinearFunction,
    NonlinearLeastSquaresProblem, NonlinearProblem, ReturnCode, isinplace,
    OverrideInit
import SciMLBase: alg_order, isadaptive
using NonlinearSolveBase: NonlinearSolveBase
using NonlinearSolveFirstOrder: NewtonRaphson
using CommonSolve: init, solve, step!
using DiffEqBase: DefaultInit
import DiffEqBase: initialize!
using OrdinaryDiffEqCore: OrdinaryDiffEqCore, OrdinaryDiffEqAlgorithm,
    OrdinaryDiffEqMutableCache, AbstractController, AbstractControllerCache
import OrdinaryDiffEqCore: alg_cache, get_fsalfirstlast, isfsal, perform_step!,
    isdiscretecache, isdiscretealg, beta2_default, beta1_default, dt_required,
    _initialize_dae!, allows_null_u0

using Reexport: Reexport, @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("cache.jl")
include("solve.jl")
include("alg_utils.jl")
include("controller.jl")

export IDSolve

end

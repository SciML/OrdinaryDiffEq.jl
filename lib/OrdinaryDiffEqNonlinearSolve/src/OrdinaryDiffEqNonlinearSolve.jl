module OrdinaryDiffEqNonlinearSolve

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SciMLBase
import DiffEqBase
import PreallocationTools
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg

import OrdinaryDiffEq: resize_nlsolver!, _initialize_dae!
import OrdinaryDiffEqDifferentiation: update_W!, isnewton

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

end

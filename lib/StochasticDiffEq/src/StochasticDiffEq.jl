module StochasticDiffEq

using Reexport: @reexport

@reexport using DiffEqBase
@reexport using StochasticDiffEqCore
@reexport using StochasticDiffEqLowOrder
@reexport using StochasticDiffEqRODE
@reexport using StochasticDiffEqHighOrder
@reexport using StochasticDiffEqMilstein
@reexport using StochasticDiffEqROCK
@reexport using StochasticDiffEqImplicit
@reexport using StochasticDiffEqWeak
@reexport using StochasticDiffEqIIF
@reexport using StochasticDiffEqLeaping
@reexport using DiffEqNoiseProcess
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

import SciMLBase
import PrecompileTools
import Preferences

# Cross-subpackage composites (SOSRI2 from HighOrder + stiff algs from Implicit).

"""
    AutoSOSRI2(alg; kwargs...)

Automatic stiffness switching between [`SOSRI2`](@ref) and a stiff SDE algorithm
`alg` (typically an implicit method). See [`AutoAlgSwitch`](@ref).
"""
AutoSOSRI2(alg; kwargs...) = AutoAlgSwitch(SOSRI2(), alg; kwargs...)

"""
    AutoSOSRA2(alg; kwargs...)

Automatic stiffness switching between [`SOSRA2`](@ref) and a stiff SDE algorithm
`alg` (typically an implicit method). See [`AutoAlgSwitch`](@ref).
"""
AutoSOSRA2(alg; kwargs...) = AutoAlgSwitch(SOSRA2(), alg; kwargs...)

include("default_sde_alg.jl")
include("precompilation.jl")

export AutoSOSRI2, AutoSOSRA2
export NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

export solve, init, solve!, step!

export checkSRIOrder, checkSRAOrder, checkRIOrder, checkRSOrder,
    checkNONOrder,
    constructSRIW1, constructSRA1,
    constructDRI1, constructRI1, constructRI3, constructRI5, constructRI6,
    constructRDI1WM, constructRDI2WM, constructRDI3WM, constructRDI4WM,
    constructRS1, constructRS2,
    constructNON, constructNON2

end # module

module StochasticDiffEq

using Reexport
@reexport using DiffEqBase

using Reexport: @reexport

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
import OrdinaryDiffEqCore: perform_step!, loopheader!, loopfooter!

# AutoSOSRI2/AutoSOSRA2 reference concrete types from multiple solver subpackages
# (SOSRI2 from HighOrder, implicit algs from Implicit), so they live here in the umbrella.
AutoSOSRI2(alg; kwargs...) = AutoAlgSwitch(SOSRI2(), alg; kwargs...)
AutoSOSRA2(alg; kwargs...) = AutoAlgSwitch(SOSRA2(), alg; kwargs...)

include("default_sde_alg.jl")

export AutoSOSRI2, AutoSOSRA2
export NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

# Re-export general functions
export solve, init, solve!, step!

# Misc Tools (tableau constructors from HighOrder)
export checkSRIOrder, checkSRAOrder, checkRIOrder, checkRSOrder,
    checkNONOrder,
    constructSRIW1, constructSRA1,
    constructDRI1, constructRI1, constructRI3, constructRI5, constructRI6,
    constructRDI1WM, constructRDI2WM, constructRDI3WM, constructRDI4WM,
    constructRS1, constructRS2,
    constructNON, constructNON2

end # module

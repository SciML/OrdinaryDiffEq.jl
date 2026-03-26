"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

import DocStringExtensions
using Reexport

# OrdinaryDiffEqCore provides core utilities and re-exports SciMLBase.
@reexport using OrdinaryDiffEqCore

# OrdinaryDiffEqDefault provides: DefaultODEAlgorithm and the default algorithm
# selection logic. It transitively depends on Tsit5, Verner, Rosenbrock, and BDF.
@reexport using OrdinaryDiffEqDefault

# Re-export the commonly used algorithms from OrdinaryDiffEqDefault's dependencies.
# These are the algorithms that the default algorithm selector can choose from.
@reexport using OrdinaryDiffEqTsit5
@reexport using OrdinaryDiffEqVerner
@reexport using OrdinaryDiffEqRosenbrock
@reexport using OrdinaryDiffEqBDF

# Re-export specific core utilities that were previously exported
using OrdinaryDiffEqCore: OrdinaryDiffEqCore,
    CompositeAlgorithm,
    ShampineCollocationInit, BrownFullBasicInit, NoInit,
    du_cache, full_cache, isfsal, ode_interpolant, u_cache,
    AutoSwitch,
    IController, PIController, PIDController

import CommonSolve: init, solve, solve!, step!
import SciMLBase: SciMLBase, addsteps!, savevalues!, terminate!

# General Functions
export solve, solve!, init, step!

# Callback Necessary
export addsteps!, ode_interpolant, terminate!, savevalues!, isfsal

# Core types
export CompositeAlgorithm, AutoSwitch
export ShampineCollocationInit, BrownFullBasicInit, NoInit
export IController, PIController, PIDController

# Re-export Reexport for downstream compatibility
export Reexport, @reexport

end # module

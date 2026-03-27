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

# Import only the widely-used algorithms from OrdinaryDiffEqDefault's dependencies.
# We intentionally do NOT @reexport the full sub-packages so that less common
# algorithms (e.g. other Rosenbrock or BDF variants) stay opt-in.
using OrdinaryDiffEqTsit5: Tsit5, AutoTsit5
using OrdinaryDiffEqVerner: Vern6, Vern7, Vern8, Vern9,
    AutoVern6, AutoVern7, AutoVern8, AutoVern9
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
using OrdinaryDiffEqBDF: FBDF

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

# Widely-used algorithms (selectively exported, not blanket @reexport)
export Tsit5, AutoTsit5
export Vern6, Vern7, Vern8, Vern9, AutoVern6, AutoVern7, AutoVern8, AutoVern9
export Rosenbrock23, Rodas5P
export FBDF

# Re-export Reexport for downstream compatibility
export Reexport, @reexport

end # module

module OrdinaryDiffEqSymplecticRK

import OrdinaryDiffEqCore: perform_step!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqPartitionedAlgorithm,
    CompiledFloats,
    alg_cache, @cache,
    constvalue,
    get_fsalfirstlast,
    generic_solver_docstring, default_linear_interpolation
import SciMLBase: alg_order, full_cache
import DiffEqBase: initialize!
using FastBroadcast: @..
using MuladdMacro: @muladd
using RecursiveArrayTools: ArrayPartition
import OrdinaryDiffEqCore

using Reexport: Reexport, @reexport
using SciMLBase: SciMLBase
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("symplectic_caches.jl")
include("symplectic_tableaus.jl")
include("symplectic_perform_step.jl")

export SymplecticEuler, VelocityVerlet, VerletLeapfrog, LeapfrogDriftKickDrift,
    PseudoVerletLeapfrog, McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
    CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

end

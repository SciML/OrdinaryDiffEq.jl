module OrdinaryDiffEqSymplecticRK

import OrdinaryDiffEq: OrdinaryDiffEqPartitionedAlgorithm, OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       @cache, @unpack
import MuladdMacro: @muladd


include("algorithms.jl")
include("alg_utils.jl")
include("symplectic_caches.jl")
include("symplectic_perform_step.jl")
include("symplectic_tableaus.jl")

export SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog,
       McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
       CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

end

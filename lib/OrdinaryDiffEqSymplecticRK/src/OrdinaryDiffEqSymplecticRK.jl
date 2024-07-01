module OrdinaryDiffEqSymplecticRK

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant!, _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools


include("algorithms.jl")
include("alg_utils.jl")
include("symplectic_caches.jl")
include("symplectic_perform_step.jl")
include("symplectic_tableaus.jl")

export SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog,
       McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
       CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

end

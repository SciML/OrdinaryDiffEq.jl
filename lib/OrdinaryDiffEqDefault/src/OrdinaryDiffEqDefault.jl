module OrdinaryDiffEqDefault

using OrdinaryDiffEqCore: alg_stability_size, beta2_default, beta1_default, AutoSwitchCache,
                      ODEIntegrator, trivial_limiter!,
                      CompositeAlgorithm, OrdinaryDiffEqAlgorithm,
                      OrdinaryDiffEqMutableCache, AutoAlgSwitch
using OrdinaryDiffEqVerner: Vern7, Vern8, Vern9, Vern6
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
using OrdinaryDiffEqBDF: FBDF

import OrdinaryDiffEqCore: is_mass_matrix_alg, default_autoswitch
import LinearSolve
using LinearAlgebra: I, isdiag
using EnumX

using Reexport
@reexport using DiffEqBase

include("default_alg.jl")

export DefaultODEAlgorithm

end # module OrdinaryDiffEqDefault

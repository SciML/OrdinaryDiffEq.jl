module OrdinaryDiffEqDefault

using OrdinaryDiffEq: Vern7, Vern8, Vern9, Vern6, Tsit5, FBDF,
    alg_stability_size, beta2_default, beta1_default, AutoSwitchCache, ODEIntegrator,
    CompositeAlgorithm, OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache, AutoAlgSwitch
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
import OrdinaryDiffEq: is_mass_matrix_alg, default_autoswitch
import LinearSolve
using LinearAlgebra: I, isdiag
using EnumX

include("default_alg.jl")

export DefaultODEAlgorithm

end # module OrdinaryDiffEqDefault

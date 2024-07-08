module OrdinaryDiffEqDefault

using OrdinaryDiffEq: Vern7, Vern8, Vern9, Vern6, Tsit5, Rosenbrock23, Rodas5P, FBDF,
    alg_stability_size, beta2_default, beta1_default, AutoSwitchCache, ODEIntegrator,
    CompositeAlgorithm, OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache, AutoAlgSwitch
import OrdinaryDiffEq: is_mass_matrix_alg, default_autoswitch
import LinearSolve
using LinearAlgebra: I
using EnumX

include("default_alg.jl")

export DefaultODEAlgorithm

end # module OrdinaryDiffEqDefault

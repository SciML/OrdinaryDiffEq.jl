module OrdinaryDiffEqAMF

import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm

using LinearAlgebra
using LinearSolve
using SciMLOperators

using Reexport
@reexport using SciMLBase

include("linsolve.jl")
include("operators.jl")
include("amf.jl")

export AMF
export SciMLOpFactorization
export AMFOperator
export split_jacobian_operator
export build_amf_function

end # module OrdinaryDiffEqAMF

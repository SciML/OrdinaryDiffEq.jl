module OrdinaryDiffEqFeagin

include("alg_utils.jl")
include("algorithms.jl")
include("feagin_caches.jl")
include("feagin_rk_perform_step.jl")
include("feagin_tableaus.jl")

export Feagin10, Feagin12, Feagin14

end
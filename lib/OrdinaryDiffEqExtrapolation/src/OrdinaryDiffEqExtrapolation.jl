module OrdinaryDiffEqExtrapolation

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("extrapolation_caches.jl")
include("extrapolation_perform_step.jl")

@inline function DiffEqBase.get_tmp_cache(integrator,
    alg::OrdinaryDiffEqImplicitExtrapolationAlgorithm,
    cache::OrdinaryDiffEqMutableCache)
(cache.tmp, cache.utilde)
end

export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
       ImplicitEulerExtrapolation,
       ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
       ImplicitEulerBarycentricExtrapolation

end
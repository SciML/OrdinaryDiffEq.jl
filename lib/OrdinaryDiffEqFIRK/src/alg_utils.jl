qmax_default(alg::Union{RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau}) = 8

alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::RadauIIA9) = 9
alg_order(alg::AdaptiveRadau) = 5 #dummy value

isfirk(alg::RadauIIA3) = true
isfirk(alg::RadauIIA5) = true
isfirk(alg::RadauIIA9) = true
isfirk(alg::AdaptiveRadau) = true

alg_adaptive_order(alg::RadauIIA3) = 1
alg_adaptive_order(alg::RadauIIA5) = 3
alg_adaptive_order(alg::RadauIIA9) = 5

get_current_alg_order(alg::AdaptiveRadau, cache) = cache.num_stages * 2 - 1
get_current_adaptive_order(alg::AdaptiveRadau, cache) = cache.num_stages

function has_stiff_interpolation(::Union{RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau})
    return true
end

# Type-stable controller trait dispatches
for Alg in [:RadauIIA3, :RadauIIA5, :RadauIIA9, :AdaptiveRadau]
    @eval ispredictive(alg::$Alg) = ispredictive(alg.controller)
    @eval isstandard(alg::$Alg) = isstandard(alg.controller)
end

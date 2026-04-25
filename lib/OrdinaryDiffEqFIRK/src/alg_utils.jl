qmax_default(alg::Union{RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau}) = 8

alg_order(alg::RadauIIA3) = 3
alg_order(alg::RadauIIA5) = 5
alg_order(alg::RadauIIA9) = 9
alg_order(alg::AdaptiveRadau) = 5 #dummy value

default_controller(QT, alg::RadauIIA3) = PredictiveController(QT, alg)
default_controller(QT, alg::RadauIIA5) = PredictiveController(QT, alg)
default_controller(QT, alg::RadauIIA9) = PredictiveController(QT, alg)
default_controller(QT, alg::AdaptiveRadau) = PredictiveController(QT, alg)

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

qmax_default(alg::GaussLegendre) = 8

alg_order(alg::GaussLegendre) = 2 * alg.num_stages

default_controller(QT, alg::GaussLegendre) = PIController(QT, alg)

isfirk(alg::GaussLegendre) = true

# Richardson step-doubling controller
isadaptive(alg::GaussLegendre) = alg.num_stages >= 2
alg_adaptive_order(alg::GaussLegendre) = 2 * alg.num_stages
has_stiff_interpolation(::GaussLegendre) = true

get_current_alg_order(alg::GaussLegendre, cache) = 2 * alg.num_stages
get_current_adaptive_order(alg::GaussLegendre, cache) = 2 * alg.num_stages

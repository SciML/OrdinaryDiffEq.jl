alg_order(alg::AitkenNeville) = alg.init_order
alg_maximum_order(alg::ExtrapolationMidpointDeuflhard) = 2(alg.max_order + 1)

function get_current_adaptive_order(
        alg::OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm,
        cache
    )
    return cache.cur_order
end
function get_current_alg_order(
        alg::OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm,
        cache
    )
    return cache.cur_order
end
get_current_alg_order(alg::ExtrapolationMidpointDeuflhard, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitDeuflhardExtrapolation, cache) = 2(cache.n_curr + 1)
get_current_adaptive_order(alg::ExtrapolationMidpointDeuflhard, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitDeuflhardExtrapolation, cache) = 2cache.n_curr
get_current_alg_order(alg::ExtrapolationMidpointHairerWanner, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitHairerWannerExtrapolation, cache) = 2(cache.n_curr + 1)
get_current_alg_order(alg::ImplicitEulerBarycentricExtrapolation, cache) = cache.n_curr
get_current_alg_order(alg::ImplicitEulerExtrapolation, cache) = cache.n_curr + 1
get_current_adaptive_order(alg::ExtrapolationMidpointHairerWanner, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitHairerWannerExtrapolation, cache) = 2cache.n_curr
get_current_adaptive_order(alg::ImplicitEulerExtrapolation, cache) = cache.n_curr - 1
function get_current_adaptive_order(
        alg::ImplicitEulerBarycentricExtrapolation, cache
    )
    return cache.n_curr - 2
end

alg_maximum_order(alg::ImplicitDeuflhardExtrapolation) = 2(alg.max_order + 1)
alg_maximum_order(alg::ExtrapolationMidpointHairerWanner) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitHairerWannerExtrapolation) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitEulerExtrapolation) = 2(alg.max_order + 1)
alg_maximum_order(alg::ImplicitEulerBarycentricExtrapolation) = alg.max_order

function default_controller(
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
            ExtrapolationMidpointHairerWanner,
            ImplicitHairerWannerExtrapolation,
            ImplicitEulerExtrapolation,
            ImplicitEulerBarycentricExtrapolation,
        },
        cache,
        qoldinit, _beta1 = nothing, _beta2 = nothing
    )
    QT = typeof(qoldinit)
    beta1, beta2 = _digest_beta1_beta2(alg, cache, Val(QT), _beta1, _beta2)
    return ExtrapolationController(beta1)
end

beta2_default(alg::ExtrapolationMidpointDeuflhard) = 0 // 1
beta2_default(alg::ImplicitDeuflhardExtrapolation) = 0 // 1
beta2_default(alg::ExtrapolationMidpointHairerWanner) = 0 // 1
beta2_default(alg::ImplicitHairerWannerExtrapolation) = 0 // 1
beta2_default(alg::ImplicitEulerExtrapolation) = 0 // 1
beta2_default(alg::ImplicitEulerBarycentricExtrapolation) = 0 // 1

beta1_default(alg::ExtrapolationMidpointDeuflhard, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitDeuflhardExtrapolation, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ExtrapolationMidpointHairerWanner, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitHairerWannerExtrapolation, beta2) = 1 // (2alg.init_order + 1)
beta1_default(alg::ImplicitEulerExtrapolation, beta2) = 1 // (alg.init_order + 1)
beta1_default(alg::ImplicitEulerBarycentricExtrapolation, beta2) = 1 // (alg.init_order - 1)

function gamma_default(alg::ExtrapolationMidpointDeuflhard)
    return (1 // 4)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitDeuflhardExtrapolation)
    return (1 // 4)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ExtrapolationMidpointHairerWanner)
    return (65 // 100)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitHairerWannerExtrapolation)
    return (65 // 100)^beta1_default(alg, beta2_default(alg))
end
function gamma_default(alg::ImplicitEulerExtrapolation)
    return (65 // 100)^beta1_default(alg, beta2_default(alg))
end

function gamma_default(alg::ImplicitEulerBarycentricExtrapolation)
    return (80 // 100)^beta1_default(alg, beta2_default(alg))
end

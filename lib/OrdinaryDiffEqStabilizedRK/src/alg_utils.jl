alg_order(alg::ROCK2) = 2
alg_order(alg::ROCK4) = 4

alg_order(alg::ESERK4) = 4
alg_order(alg::ESERK5) = 5
alg_order(alg::SERK2) = 2

alg_order(alg::RKC) = 2

alg_order(alg::TSRKC3) = 3

alg_extrapolates(alg::TSRKC3) = true

default_controller(QT, alg::SERK2) = PredictiveController(QT, alg)
default_controller(QT, alg::Union{RKC, TSRKC3}) = PredictiveController(QT, alg)
default_controller(QT, alg::Union{RKL1, RKL2}) = PredictiveController(QT, alg)

alg_adaptive_order(alg::RKL1) = 1
alg_adaptive_order(alg::RKL2) = 2
alg_adaptive_order(alg::Union{RKC, TSRKC3}) = 2

gamma_default(alg::Union{RKC, TSRKC3}) = 8 // 10
gamma_default(alg::Union{RKL1, RKL2}) = 8 // 10
qmax_default(alg::TSRKC3) = 2

fac_default_gamma(alg::Union{RKC, SERK2, TSRKC3}) = true
fac_default_gamma(alg::Union{RKL1, RKL2}) = true
has_dtnew_modification(alg::Union{ROCK2, ROCK4, SERK2, ESERK4, ESERK5}) = true

function dtnew_modification(integrator, alg::ROCK2, dtnew)
    return min(
        dtnew,
        typeof(dtnew)(
            (
                (
                    (min(integrator.alg.max_stages, 200)^2) * 0.811 -
                        1.5
                ) / integrator.eigen_est
            )
        )
    )
end
function dtnew_modification(integrator, alg::ROCK4, dtnew)
    return min(
        dtnew,
        typeof(dtnew)(
            (
                ((min(integrator.alg.max_stages, 152)^2) * 0.353 - 3) /
                    integrator.eigen_est
            )
        )
    )
end
function dtnew_modification(integrator, alg::SERK2, dtnew)
    return min(
        dtnew,
        typeof(dtnew)((0.8 * 250 * 250 / (integrator.eigen_est + 1)))
    )
end
function dtnew_modification(integrator, alg::ESERK4, dtnew)
    return min(dtnew, typeof(dtnew)((0.98 * 4000 * 4000 / integrator.eigen_est)))
end
function dtnew_modification(integrator, alg::ESERK5, dtnew)
    return min(dtnew, typeof(dtnew)((0.98 * 2000 * 2000 / integrator.eigen_est)))
end


alg_order(alg::RKL1) = 1
alg_order(alg::RKL2) = 2

# clamp stage counts to odd ints that are plausible
@inline function _rkl_clamp_odd_stages(min_stages::Int, max_stages::Int)
    min_stage = max(3, min_stages)
    min_stage = isodd(min_stage) ? min_stage : min_stage + 1
    max_stage = max(max_stages, min_stage)
    max_stage = isodd(max_stage) ? max_stage : max_stage - 1
    if max_stage < min_stage
        max_stage = min_stage
    end
    return min_stage, max_stage
end

# cap dt so that even at max stage count the stability bound is not exceeded
has_dtnew_modification(alg::Union{RKL1, RKL2}) = true

function dtnew_modification(integrator, alg::RKL1, dtnew)
    _, s = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    return min(dtnew, typeof(dtnew)((s^2 + s) / integrator.eigen_est))
end

function dtnew_modification(integrator, alg::RKL2, dtnew)
    _, s = _rkl_clamp_odd_stages(alg.min_stages, alg.max_stages)
    return min(dtnew, typeof(dtnew)((s^2 + s - 2) / (2 * integrator.eigen_est)))
end

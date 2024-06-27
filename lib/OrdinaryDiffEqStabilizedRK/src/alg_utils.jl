alg_order(alg::ROCK2) = 2
alg_order(alg::ROCK4) = 4

alg_order(alg::ESERK4) = 4
alg_order(alg::ESERK5) = 5
alg_order(alg::SERK2) = 2

alg_order(alg::RKC) = 2

ispredictive(alg::Union{SERK2}) = alg.controller === :Predictive
ispredictive(alg::Union{RKC}) = true

alg_adaptive_order(alg::RKC) = 2

gamma_default(alg::RKC) = 8 // 10

fac_default_gamma(alg::Union{RKC, SERK2}) = true
has_dtnew_modification(alg::Union{ROCK2, ROCK4, SERK2, ESERK4, ESERK5}) = true

function dtnew_modification(integrator, alg::ROCK2, dtnew)
    min(dtnew,
        typeof(dtnew)((((min(integrator.alg.max_stages, 200)^2.0) * 0.811 -
                        1.5) / integrator.eigen_est)))
end
function dtnew_modification(integrator, alg::ROCK4, dtnew)
    min(dtnew,
        typeof(dtnew)((((min(integrator.alg.max_stages, 152)^2.0) * 0.353 - 3) /
                       integrator.eigen_est)))
end
function dtnew_modification(integrator, alg::SERK2, dtnew)
    min(dtnew,
        typeof(dtnew)((0.8 * 250 * 250 / (integrator.eigen_est + 1.0))))
end
function dtnew_modification(integrator, alg::ESERK4, dtnew)
    min(dtnew, typeof(dtnew)((0.98 * 4000 * 4000 / integrator.eigen_est)))
end
function dtnew_modification(integrator, alg::ESERK5, dtnew)
    min(dtnew, typeof(dtnew)((0.98 * 2000 * 2000 / integrator.eigen_est)))
end

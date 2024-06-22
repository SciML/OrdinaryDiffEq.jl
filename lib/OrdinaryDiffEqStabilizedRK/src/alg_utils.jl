alg_order(alg::ROCK2) = 2
alg_order(alg::ROCK4) = 4

alg_order(alg::ESERK4) = 4
alg_order(alg::ESERK5) = 5
alg_order(alg::SERK2) = 2

alg_order(alg::IRKC) = 2
alg_order(alg::RKC) = 2

ispredictive(alg::Union{SERK2}) = alg.controller === :Predictive
ispredictive(alg::Union{RKC}) = true

alg_adaptive_order(alg::RKC) = 2
alg_adaptive_order(alg::IRKC) = 1

gamma_default(alg::RKC) = 8 // 10
gamma_default(alg::IRKC) = 8 // 10

alg_can_repeat_jac(alg::IRKC) = false

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
alg_order(alg::RandomEM) = 1 // 2
alg_order(alg::RandomHeun) = 1 // 2
alg_order(alg::RandomTamedEM) = 1 // 2
alg_order(alg::BAOAB) = 1 // 1

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BAOAB) = is_diagonal_noise(prob)

alg_order(alg::TanYam7) = 7
alg_order(alg::DP8) = 8
alg_order(alg::TsitPap8) = 8
alg_order(alg::PFRK87) = 8

qmin_default(alg::DP8) = 1 // 3

qmax_default(alg::DP8) = 6

beta2_default(alg::DP8) = 0 // 1

beta1_default(alg::DP8, beta2) = typeof(beta2)(1 // alg_order(alg)) - beta2 / 5
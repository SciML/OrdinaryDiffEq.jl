alg_order(alg::AB3) = 3
alg_order(alg::ABM32) = 3
alg_order(alg::AB4) = 4
alg_order(alg::ABM43) = 4
alg_order(alg::AB5) = 5
alg_order(alg::ABM54) = 5
alg_order(alg::VCAB3) = 3
alg_order(alg::VCAB4) = 4
alg_order(alg::VCAB5) = 5
alg_order(alg::VCABM3) = 3
alg_order(alg::VCABM4) = 4
alg_order(alg::VCABM5) = 5
alg_order(alg::VCABM) = 1  #dummy value

isstandard(alg::VCABM) = true
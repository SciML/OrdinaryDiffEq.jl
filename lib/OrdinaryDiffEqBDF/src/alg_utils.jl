alg_extrapolates(alg::ABDF2) = true
alg_extrapolates(alg::SBDF) = true
alg_extrapolates(alg::MEBDF2) = true

alg_order(alg::ABDF2) = 2
alg_order(alg::SBDF) = alg.order
alg_order(alg::QNDF1) = 1
alg_order(alg::MEBDF2) = 2
alg_order(alg::FBDF) = 1 #dummy value
alg_order(alg::QNDF2) = 2
alg_order(alg::QNDF) = 1 #dummy value

get_current_alg_order(alg::QNDF, cache) = cache.order
get_current_alg_order(alg::FBDF, cache) = cache.order

get_current_adaptive_order(alg::FBDF, cache) = cache.order

qsteady_max_default(alg::QNDF) = 2 // 1
qsteady_max_default(alg::QNDF2) = 2 // 1
qsteady_max_default(alg::QNDF1) = 2 // 1

qsteady_min_default(alg::FBDF) = 9 // 10
qsteady_max_default(alg::FBDF) = 2 // 1
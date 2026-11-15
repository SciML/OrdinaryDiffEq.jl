alg_order(alg::IRKC) = 2
gamma_default(alg::IRKC) = 8 // 10
alg_can_repeat_jac(alg::IRKC) = false
issplit(alg::IRKC) = true
fac_default_gamma(alg::IRKC) = true
default_controller(QT, alg::IRKC) = IController(QT, alg)

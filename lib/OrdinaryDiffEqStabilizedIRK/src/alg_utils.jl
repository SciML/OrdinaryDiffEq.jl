alg_order(alg::IRKC) = 2
alg_adaptive_order(alg::IRKC) = 1
gamma_default(alg::IRKC) = 8 // 10
alg_can_repeat_jac(alg::IRKC) = false
issplit(alg::IRKC) = true
fac_default_gamma(alg::IRKC) = true

# Type-stable default_controller_v7 dispatch for IRKC
@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    default_controller_v7(QT, alg::IRKC) = NewIController(QT, alg)
end

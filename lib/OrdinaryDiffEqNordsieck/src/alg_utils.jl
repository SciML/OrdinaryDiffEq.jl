alg_order(alg::AN5) = 5
alg_order(alg::JVODE) = 1  #dummy value

alg_adaptive_order(alg::AN5) = 5
qsteady_max_default(alg::AN5) = 3 // 2

# JVODE-tuned CommonControllerOptions defaults match the historical alg-struct kwargs
# `JVODE(qmin=1/5, qmax=10, qsteady_min=1, qsteady_max=3//2)`.
qmin_default(::JVODE) = 1 // 5
qmax_default(::JVODE) = 10 // 1
qsteady_min_default(::JVODE) = 1 // 1
qsteady_max_default(::JVODE) = 3 // 2

get_current_alg_order(alg::JVODE, cache) = get_current_adaptive_order(alg, cache)

function OrdinaryDiffEqCore.default_controller(QT, alg::JVODE)
    # Thread the alg-struct kwargs through to the controller so that
    # `JVODE(qmax = 20)` keeps working.
    return JVODEController(
        QT, alg;
        qmin = alg.qmin, qmax = alg.qmax,
        qsteady_min = alg.qsteady_min, qsteady_max = alg.qsteady_max,
    )
end

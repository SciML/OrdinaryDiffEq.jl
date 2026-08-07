alg_extrapolates(alg::ABDF2) = true
alg_extrapolates(alg::SBDF) = true
alg_extrapolates(alg::MEBDF2) = true

alg_order(alg::ABDF2) = 2
alg_order(alg::SBDF) = alg.order
alg_order(alg::QNDF1) = 1
alg_order(alg::QNDF2) = 2
alg_order(alg::QNDF) = 1 #dummy value
alg_order(alg::MEBDF2) = 2
alg_order(alg::FBDF) = 1 #dummy value

issplit(alg::SBDF) = true

qsteady_min_default(alg::FBDF) = 9 // 10

qsteady_max_default(alg::QNDF) = 2 // 1
qsteady_max_default(alg::QNDF2) = 2 // 1
qsteady_max_default(alg::QNDF1) = 2 // 1
qsteady_max_default(alg::FBDF) = 2 // 1

get_current_alg_order(alg::QNDF, cache) = cache.order
get_current_alg_order(alg::FBDF, cache) = cache.order

get_current_adaptive_order(alg::QNDF, cache) = cache.order
get_current_adaptive_order(alg::FBDF, cache) = cache.order

isadaptive(alg::DImplicitEuler) = true
isadaptive(alg::DABDF2) = true
isadaptive(alg::DFBDF) = true

has_special_newton_error(alg::QNDF) = true

default_controller(QT, alg::DImplicitEuler) = IController(QT, alg)
default_controller(QT, alg::ABDF2) = IController(QT, alg)
default_controller(QT, alg::QNDF1) = IController(QT, alg)
default_controller(QT, alg::QNDF2) = IController(QT, alg)
default_controller(QT, alg::DABDF2) = IController(QT, alg)

alg_extrapolates(alg::DImplicitEuler) = true
alg_extrapolates(alg::DABDF2) = true

alg_order(alg::DImplicitEuler) = 1
alg_order(alg::DABDF2) = 2
alg_order(alg::DFBDF) = 1 #dummy value

isfsal(alg::DImplicitEuler) = false

has_stiff_interpolation(::Union{QNDF, FBDF, DFBDF}) = true

############################################ NordsieckBDF / DNordsieckBDF
alg_order(alg::NordsieckBDF) = 1  # dummy: the running order lives in the cache
alg_order(alg::DNordsieckBDF) = 1
isadaptive(alg::DNordsieckBDF) = true
get_current_alg_order(alg::NordsieckBDFAlgs, cache) = cache.order
get_current_adaptive_order(alg::NordsieckBDFAlgs, cache) = cache.order
has_stiff_interpolation(::NordsieckBDFAlgs) = true

# The Newton increment norm is scaled by tq[2], which converts it into the units of
# the local error test. That makes `NLNewton(κ = ...)` mean CVODE's NLSCOEF: the
# fraction of the error-test budget the corrector may consume.
has_special_newton_error(alg::NordsieckBDFAlgs) = true


# The step-size logic is CVODE's (`cvSetEta` keeps h unless eta >= 1.5), so the

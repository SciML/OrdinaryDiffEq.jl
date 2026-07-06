alg_order(alg::TauLeaping) = 1 // 1
alg_order(alg::CaoTauLeaping) = 1 // 1
alg_order(alg::ImplicitTauLeaping) = 1 // 1  # Weak first order (backward Euler)
alg_order(alg::ThetaTrapezoidalTauLeaping) = 2 // 1  # Weak second order

# special cases in stepsize_controllers.jl
# Match OrdinaryDiffEqCore's `default_controller(QT, alg)` argument order (QT first,
# alg second). The previous `(alg, args...)` slurp put `alg` in the `QT` position, so
# `solve`'s `default_controller(QT, alg)` never reached it (it fell through to the
# generic `PIController` method) and it was ambiguous with
# `default_controller(QT, ::OrdinaryDiffEqCompositeAlgorithm)`.
function OrdinaryDiffEqCore.default_controller(
        QT, alg::Union{TauLeaping, CaoTauLeaping, ThetaTrapezoidalTauLeaping}
    )
    return DummyController()
end

# CaoTauLeaping always accepts (a posteriori dt), others use error estimates.
OrdinaryDiffEqCore.isaposteriori(alg::CaoTauLeaping) = true

alg_control_rate(::TauLeaping) = true
alg_control_rate(::ThetaTrapezoidalTauLeaping) = true  # We manage rates manually

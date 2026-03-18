alg_order(alg::TauLeaping) = 1 // 1
alg_order(alg::CaoTauLeaping) = 1 // 1
alg_order(alg::ImplicitTauLeaping) = 1 // 1  # Weak first order (backward Euler)
alg_order(alg::ThetaTrapezoidalTauLeaping) = 2 // 1  # Weak second order

# special cases in stepsize_controllers.jl
function StochasticDiffEqCore.default_controller(
        alg::Union{TauLeaping, CaoTauLeaping, ThetaTrapezoidalTauLeaping}, args...
    )
    return DummyController()
end

# CaoTauLeaping always accepts (a posteriori dt), others use error estimates.
OrdinaryDiffEqCore.isaposteriori(alg::CaoTauLeaping) = true

alg_control_rate(::TauLeaping) = true
alg_control_rate(::ThetaTrapezoidalTauLeaping) = true  # We manage rates manually

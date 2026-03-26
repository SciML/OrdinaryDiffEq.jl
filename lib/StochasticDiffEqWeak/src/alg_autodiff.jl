# _alg_autodiff methods for SDE Newton algorithm abstract types.
# These extract the autodiff field from the algorithm struct so that
# OrdinaryDiffEqDifferentiation can set up the Jacobian computation correctly.

OrdinaryDiffEqDifferentiation._alg_autodiff(alg::StochasticDiffEqNewtonAlgorithm) = alg.autodiff
OrdinaryDiffEqDifferentiation._alg_autodiff(alg::StochasticDiffEqNewtonAdaptiveAlgorithm) = alg.autodiff
OrdinaryDiffEqDifferentiation._alg_autodiff(alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm) = alg.autodiff
OrdinaryDiffEqDifferentiation._alg_autodiff(alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm) = alg.autodiff

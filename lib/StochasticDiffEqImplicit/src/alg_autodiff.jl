# _alg_autodiff methods for SDE Newton algorithm abstract types.
# These extract the AD type parameter from the algorithm struct so that
# OrdinaryDiffEqDifferentiation can set up the Jacobian computation correctly.

function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end

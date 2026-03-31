struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split = true) = EM{split}()

struct SplitEM <: StochasticDiffEqAlgorithm end

struct EulerHeun <: StochasticDiffEqAlgorithm end

struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split = true) = LambaEM{split}()

struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

struct SimplifiedEM <: StochasticDiffEqAlgorithm end

struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(; interpretation = SciMLBase.AlgorithmInterpretation.Ito) = RKMil{interpretation}()

struct RKMilCommute{T} <: StochasticDiffEqAdaptiveAlgorithm
    interpretation::SciMLBase.AlgorithmInterpretation.T
    ii_approx::T
end
function RKMilCommute(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, ii_approx = IICommutative())
    return RKMilCommute(interpretation, ii_approx)
end

struct PCEuler{T <: Real, F} <: StochasticDiffEqAlgorithm
    theta::T
    eta::T
    ggprime::F
end
PCEuler(ggprime; theta = 1 / 2, eta = 1 / 2) = PCEuler(theta, eta, ggprime)

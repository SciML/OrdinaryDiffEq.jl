"""
@article{feagin2012high,
title={High-order explicit Runge-Kutta methods using m-symmetry},
author={Feagin, Terry},
year={2012},
publisher={Neural, Parallel \\& Scientific Computations}
}

Feagin10: Explicit Runge-Kutta Method
Feagin's 10th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin10{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
@article{feagin2012high,
title={High-order explicit Runge-Kutta methods using m-symmetry},
author={Feagin, Terry},
year={2012},
publisher={Neural, Parallel \\& Scientific Computations}
}

Feagin12: Explicit Runge-Kutta Method
Feagin's 12th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin12{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
Feagin, T., “An Explicit Runge-Kutta Method of Order Fourteen,” Numerical
Algorithms, 2009

Feagin14: Explicit Runge-Kutta Method
Feagin's 14th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin14{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end
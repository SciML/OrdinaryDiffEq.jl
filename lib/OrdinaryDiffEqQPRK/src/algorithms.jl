@doc explicit_rk_docstring(
    "Runge–Kutta pairs of orders 9(8) for use in quadruple precision computations", "QPRK98",
    references = "Kovalnogov VN, Fedorov RV, Karpukhina TV, Simos TE, Tsitouras C. Runge–Kutta pairs 
    of orders 9 (8) for use in quadruple precision computations. Numerical Algorithms, 2023. 
    doi: https://doi.org/10.1007/s11075-023-01632-8")
Base.@kwdef struct QPRK98{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function QPRK98(stage_limiter!, step_limiter! = trivial_limiter!)
    QPRK98(stage_limiter!, step_limiter!, False())
end
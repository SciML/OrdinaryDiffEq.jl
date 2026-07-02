@doc explicit_rk_docstring(
    "Runge–Kutta pairs of orders 9(8) for use in quadruple precision computations", "QPRK98",
    references = "Kovalnogov VN, Fedorov RV, Karpukhina TV, Simos TE, Tsitouras C. Runge–Kutta pairs 
    of orders 9 (8) for use in quadruple precision computations. Numerical Algorithms, 2023. 
    doi: https://doi.org/10.1007/s11075-023-01632-8"
)
Base.@kwdef struct QPRK98{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
    # When true, the cache allocates dedicated rate buffers so `init`/`solve`/
    # `reinit!` compute the starting `dt` without allocating. Default false:
    # the cache footprint is unchanged from the historical layout and `initdt`
    # allocates its rate temporaries at call time (the state/unit-less scratch
    # is reused either way).
    preallocate_initdt_buffers::Bool = false
end

OrdinaryDiffEqCore.has_stage_limiter(::QPRK98) = true

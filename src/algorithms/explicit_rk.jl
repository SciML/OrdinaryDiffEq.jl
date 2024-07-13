function Base.show(io::IO, alg::OrdinaryDiffEqAlgorithm)
    print(io, String(typeof(alg).name.name), "(;")
    for fieldname in fieldnames(typeof(alg))
        print(io, " ", fieldname, " = ", getfield(alg, fieldname), ",")
    end
    print(io, ")")
end

"""
Euler - The canonical forward Euler method. Fixed timestep only.
"""

"""
KuttaPRK2p5: Parallel Explicit Runge-Kutta Method
A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.

These methods utilize multithreading on the f calls to parallelize the problem.
This requires that simultaneous calls to f are thread-safe.
"""
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end

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

"""
KuttaPRK2p5: Parallel Explicit Runge-Kutta Method
A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.

These methods utilize multithreading on the f calls to parallelize the problem.
This requires that simultaneous calls to f are thread-safe.
"""
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end
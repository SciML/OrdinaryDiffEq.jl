@doc explicit_rk_docstring(
    "Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant). Recommended for most non-stiff problems. Good default choice for unknown stiffness. Highly efficient and generic. Very good performance for most non-stiff ODEs. Recommended as default method for unknown stiffness problems.",
    "Tsit5",
    references = "@article{tsitouras2011runge,
    title={Runge--Kutta pairs of order 5 (4) satisfying only the first column simplifying assumption},
    author={Tsitouras, Ch},
    journal={Computers \\& Mathematics with Applications},
    volume={62},
    number={2},
    pages={770--775},
    year={2011},
    publisher={Elsevier},
    doi={10.1016/j.camwa.2011.06.002}
    }"
)
Base.@kwdef struct Tsit5{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
    # When true (default), the cache preallocates the scratch buffers the
    # initial-dt estimator needs, so `init`/`solve` compute the starting `dt`
    # without allocating. Set false for a leaner cache; `initdt` then allocates.
    precompute_initdt_cache::Bool = true
end
@truncate_stacktrace Tsit5 3

"""
Automatic switching algorithm that can switch between the (non-stiff) `Tsit5()` and `stiff_alg`.

    AutoTsit5(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Tsit5(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoTsit5(stiff_alg; kwargs...) = AutoAlgSwitch(Tsit5(), stiff_alg; kwargs...)

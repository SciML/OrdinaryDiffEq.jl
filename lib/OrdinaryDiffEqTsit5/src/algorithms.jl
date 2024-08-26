@doc explicit_rk_docstring(
    "A fifth-order explicit Runge-Kutta method with embedded error
    estimator of Tsitouras. Free 4th order interpolant.",
    "Tsit5",
    references = "@article{tsitouras2011runge,
    title={Runge--Kutta pairs of order 5 (4) satisfying only the first column simplifying assumption},
    author={Tsitouras, Ch},
    journal={Computers \\& Mathematics with Applications},
    volume={62},
    number={2},
    pages={770--775},
    year={2011},
    publisher={Elsevier}
    }")
Base.@kwdef struct Tsit5{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
TruncatedStacktraces.@truncate_stacktrace Tsit5 3
# for backwards compatibility
function Tsit5(stage_limiter!, step_limiter! = trivial_limiter!)
    Tsit5(stage_limiter!, step_limiter!, False())
end

"""
Automatic switching algorithm that can switch between the (non-stiff) `Tsit5()` and `stiff_alg`.

    AutoTsit5(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Tsit5(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoTsit5(stiff_alg; kwargs...) = AutoAlgSwitch(Tsit5(), stiff_alg; kwargs...)

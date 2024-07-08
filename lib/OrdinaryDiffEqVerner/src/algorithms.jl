@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 6/5 Runge-Kutta method. (lazy 6th order interpolant).",
    "Vern6",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern6{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern6 3
# for backwards compatibility
function Vern6(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern6(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 7/6 Runge-Kutta method. (lazy 7th order interpolant).",
    "Vern7",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern7{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern7 3
# for backwards compatibility
function Vern7(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern7(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 8/7 Runge-Kutta method. (lazy 8th order interpolant).",
    "Vern8",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern8{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern8 3
# for backwards compatibility
function Vern8(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern8(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy9th order interpolant).",
    "Vern9",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """, extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern9{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern9 3
# for backwards compatibility
function Vern9(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern9(stage_limiter!, step_limiter!, False(), lazy)
end

AutoVern6(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern6(lazy = lazy), alg; kwargs...)
AutoVern7(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern7(lazy = lazy), alg; kwargs...)
AutoVern8(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern8(lazy = lazy), alg; kwargs...)
AutoVern9(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern9(lazy = lazy), alg; kwargs...)
@doc explicit_rk_docstring(
    "Verner's most efficient 6/5 method (lazy 6th order interpolant).",
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
    extra_keyword_default = "lazy = true"
)
Base.@kwdef struct Vern6{StageLimiter, StepLimiter, Thread, L} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::L = Val{true}()
end
@truncate_stacktrace Vern6 3
# Convert Bool lazy to Val for backwards compatibility
function Vern6(sl::SL, stl::STL, th::TH, lazy::Bool) where {SL, STL, TH}
    Vern6{SL, STL, TH, Val{lazy}}(sl, stl, th, Val{lazy}())
end
# for backwards compatibility
function Vern6(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    return Vern6(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's most efficient 7/6 method (lazy 7th order interpolant). Good for problems requiring high accuracy. Slightly more computationally expensive than Tsit5. Performance best when parameter vector remains unchanged. Recommended for high-accuracy non-stiff problems.",
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
    extra_keyword_default = "lazy = true"
)
Base.@kwdef struct Vern7{StageLimiter, StepLimiter, Thread, L} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::L = Val{true}()
end
@truncate_stacktrace Vern7 3
# Convert Bool lazy to Val for backwards compatibility
function Vern7(sl::SL, stl::STL, th::TH, lazy::Bool) where {SL, STL, TH}
    Vern7{SL, STL, TH, Val{lazy}}(sl, stl, th, Val{lazy}())
end
# for backwards compatibility
function Vern7(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    return Vern7(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's most efficient 8/7 method (lazy 8th order interpolant).",
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
    extra_keyword_default = "lazy = true"
)
Base.@kwdef struct Vern8{StageLimiter, StepLimiter, Thread, L} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::L = Val{true}()
end
@truncate_stacktrace Vern8 3
# Convert Bool lazy to Val for backwards compatibility
function Vern8(sl::SL, stl::STL, th::TH, lazy::Bool) where {SL, STL, TH}
    Vern8{SL, STL, TH, Val{lazy}}(sl, stl, th, Val{lazy}())
end
# for backwards compatibility
function Vern8(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    return Vern8(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's most efficient 9/8 method (lazy 9th order interpolant).",
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
    """, extra_keyword_default = "lazy = true"
)
Base.@kwdef struct Vern9{StageLimiter, StepLimiter, Thread, L} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::L = Val{true}()
end
@truncate_stacktrace Vern9 3
# Convert Bool lazy to Val for backwards compatibility
function Vern9(sl::SL, stl::STL, th::TH, lazy::Bool) where {SL, STL, TH}
    Vern9{SL, STL, TH, Val{lazy}}(sl, stl, th, Val{lazy}())
end
# for backwards compatibility
function Vern9(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    return Vern9(stage_limiter!, step_limiter!, False(), lazy)
end

"""
Automatic switching algorithm that can switch between the (non-stiff) `Vern6()` and `stiff_alg`.

    AutoVern6(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Vern6(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoVern6(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern6(lazy = lazy), alg; kwargs...)
"""
Automatic switching algorithm that can switch between the (non-stiff) `Vern7()` and `stiff_alg`.

    AutoVern7(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Vern7(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoVern7(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern7(lazy = lazy), alg; kwargs...)
"""
Automatic switching algorithm that can switch between the (non-stiff) `Vern8()` and `stiff_alg`.

    AutoVern8(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Vern8(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoVern8(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern8(lazy = lazy), alg; kwargs...)
"""
Automatic switching algorithm that can switch between the (non-stiff) `Vern9()` and `stiff_alg`.

    AutoVern9(stiff_alg; kwargs...)

This method is equivalent to `AutoAlgSwitch(Vern9(), stiff_alg; kwargs...)`.
To gain access to stiff algorithms you might have to install additional libraries,
such as `OrdinaryDiffEqRosenbrock`.
"""
AutoVern9(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern9(lazy = lazy), alg; kwargs...)


@doc explicit_rk_docstring(
    "Verner's RKV76.IIa 7/6 method. Most efficient 10-stage conventional pair of orders 6 and 7 with interpolants.",
    "RKV76IIa",
    references = "@misc{verner2024rkv76iia,
    title={RKV76.IIa - A 'most efficient' Runge--Kutta (10:7(6)) pair},
    author={Verner, James H},
    year={2024},
    url={https://www.sfu.ca/~jverner/RKV76.IIa.Efficient.000003389335684.240711.FLOAT6040OnWeb}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
    """,
    extra_keyword_default = "lazy = true"
)
Base.@kwdef struct RKV76IIa{StageLimiter, StepLimiter, Thread, L} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::L = Val{true}()
end
@truncate_stacktrace RKV76IIa 3
# Convert Bool lazy to Val for backwards compatibility
function RKV76IIa(sl::SL, stl::STL, th::TH, lazy::Bool) where {SL, STL, TH}
    RKV76IIa{SL, STL, TH, Val{lazy}}(sl, stl, th, Val{lazy}())
end
# for backwards compatibility
function RKV76IIa(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    return RKV76IIa(stage_limiter!, step_limiter!, False(), lazy)
end

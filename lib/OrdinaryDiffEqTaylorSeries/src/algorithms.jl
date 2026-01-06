@doc explicit_rk_docstring(
    "A second-order explicit Taylor series method.",
    "ExplicitTaylor2"
)
Base.@kwdef struct ExplicitTaylor2{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
@truncate_stacktrace ExplicitTaylor2 3
# for backwards compatibility
function ExplicitTaylor2(stage_limiter!, step_limiter! = trivial_limiter!)
    return ExplicitTaylor2(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "An arbitrary-order explicit Taylor series method.",
    "ExplicitTaylor2"
)
Base.@kwdef struct ExplicitTaylor{P, StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    order::Val{P} = Val{1}()
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

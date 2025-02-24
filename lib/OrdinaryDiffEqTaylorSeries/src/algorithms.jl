
Base.@kwdef struct ExplicitTaylor2{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
TruncatedStacktraces.@truncate_stacktrace ExplicitTaylor2 3
# for backwards compatibility
function ExplicitTaylor2(stage_limiter!, step_limiter! = trivial_limiter!)
    ExplicitTaylor2(stage_limiter!, step_limiter!, False())
end

Base.@kwdef struct ExplicitTaylor{P, StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    order::Val{P} = Val{1}()
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

Base.@kwdef struct DAETS{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
TruncatedStacktraces.@truncate_stacktrace DAETS 3
# for backwards compatibility
function DAETS(stage_limiter!, step_limiter! = trivial_limiter!)
    DAETS(stage_limiter!, step_limiter!, False())
end
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

struct DAETS{CS, AD, FDT, ST, CJ, StageLimiter, StepLimiter, Thread} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    order::CS
    adaptive::AD
    fdtype::FDT
    calck::CJ
    tableau::ST
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function DAETS(; order = 5, adaptive = false, fdtype = Val(:central),
               calck = Val(true), tableau = nothing,
               stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
               thread = False())
    DAETS(order, adaptive, fdtype, calck, tableau, stage_limiter!, step_limiter!, thread)
end
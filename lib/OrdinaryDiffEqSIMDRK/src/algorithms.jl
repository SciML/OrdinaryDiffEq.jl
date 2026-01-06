struct MER5v2{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function MER5v2(;
        stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False()
    )
    return MER5v2{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(
        stage_limiter!,
        step_limiter!,
        thread
    )
end

# for backwards compatibility
function MER5v2(stage_limiter!, step_limiter! = trivial_limiter!)
    return MER5v2{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!, False()
    )
end

struct MER6v2{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function MER6v2(;
        stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False()
    )
    return MER6v2{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(
        stage_limiter!,
        step_limiter!,
        thread
    )
end

# for backwards compatibility
function MER6v2(stage_limiter!, step_limiter! = trivial_limiter!)
    return MER6v2{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!, False()
    )
end

struct RK6v4{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function RK6v4(;
        stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False()
    )
    return RK6v4{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(
        stage_limiter!,
        step_limiter!,
        thread
    )
end

# for backwards compatibility
function RK6v4(stage_limiter!, step_limiter! = trivial_limiter!)
    return RK6v4{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!, False()
    )
end

function Base.show(io::IO, alg::Union{MER5v2, MER6v2, RK6v4})
    return print(
        io, "$(nameof(typeof(alg)))(stage_limiter! = ", alg.stage_limiter!,
        ", step_limiter! = ", alg.step_limiter!,
        ", thread = ", alg.thread, ")"
    )
end

OrdinaryDiffEqCore.alg_order(alg::MER5v2) = 5
OrdinaryDiffEqCore.alg_order(alg::MER6v2) = 6
OrdinaryDiffEqCore.alg_order(alg::RK6v4) = 6

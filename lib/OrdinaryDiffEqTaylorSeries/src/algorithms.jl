@doc explicit_rk_docstring(
    "A second-order explicit Taylor series method.",
    "ExplicitTaylor2")
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

@doc explicit_rk_docstring(
    "An arbitrary-order explicit Taylor series method.",
    "ExplicitTaylor")
Base.@kwdef struct ExplicitTaylor{P, StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    order::Val{P} = Val{1}()
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

@doc explicit_rk_docstring(
    "An adaptive order explicit Taylor series method.",
    "ExplicitTaylorAdaptiveOrder")
Base.@kwdef struct ExplicitTaylorAdaptiveOrder{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    min_order::Int = 1
    max_order::Int = 10
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

@doc explicit_rk_docstring(
    "An implicit Taylor series method.",
    "ImplicitTaylor")
struct ImplicitTaylor2{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function ImplicitTaylor2(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ImplicitTaylor2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve, precs, extrapolant, controller, step_limiter!, AD_choice)
end

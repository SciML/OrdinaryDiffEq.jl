"""
    ExplicitTaylor2(; stage_limiter! = trivial_limiter!,
        step_limiter! = trivial_limiter!, thread = Serial())

Second-order explicit Taylor series method. `ExplicitTaylor2` evaluates the first two
time derivatives of the solution with TaylorDiff.jl and advances with a fixed step size.
Supply `dt` when calling `solve`.

The right-hand side must accept TaylorDiff numbers. In particular, avoid branches or
external functions that discard or reject Taylor coefficients.

# Fields

  - `stage_limiter!`: limiter configuration retained by the common explicit-solver
    interface. The current Taylor step does not invoke this callback.
  - `step_limiter!`: deprecated algorithm-level fallback for the accepted-step limiter.
    When no solve-level `step_limiter` is supplied, a nontrivial fallback is invoked
    centrally once per accepted step.
  - `thread`: threading configuration retained by the common explicit-solver interface.
    The second-order implementation currently performs its internal broadcasts serially.

# Keywords

  - `stage_limiter! = trivial_limiter!`: stage-limiter callback with signature
    `limiter!(u, integrator, p, t)`.
  - `step_limiter! = trivial_limiter!`: deprecated algorithm-level fallback with
    signature `limiter!(u, integrator, p, t)`. When no solve-level `step_limiter` is
    supplied, a nontrivial fallback is invoked centrally once per accepted step.
  - `thread = Serial()`: FastBroadcast threading mode stored with the algorithm.

# Returns

An `ExplicitTaylor2` algorithm object for use with `solve` or `init`.

# Examples

```julia
using OrdinaryDiffEqTaylorSeries: ExplicitTaylor2
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0))
sol = solve(prob, ExplicitTaylor2(), dt = 0.05)
sol(1.0)
```
"""
Base.@kwdef struct ExplicitTaylor2{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end
@truncate_stacktrace ExplicitTaylor2 3

"""
    ExplicitTaylor(; order = Val{1}(), stage_limiter! = trivial_limiter!,
        step_limiter! = trivial_limiter!, thread = Serial())

Fixed-order explicit Taylor series method with adaptive step-size control. The method
symbolically constructs a Taylor jet of the requested order and reuses it while solving.
Choose `order = Val(N)` to use an order-`N` expansion.

The right-hand side must be traceable by Symbolics.jl and must accept TaylorDiff numbers.
Changing the right-hand-side function, parameters, or order constructs a different cached
jet, so these should remain stable across repeated solves when compilation cost matters.

# Fields

  - `order`: compile-time `Val` containing the Taylor expansion order.
  - `stage_limiter!`: limiter configuration retained by the common explicit-solver
    interface. The current Taylor step does not invoke this callback.
  - `step_limiter!`: deprecated algorithm-level fallback for the accepted-step limiter.
    When no solve-level `step_limiter` is supplied, a nontrivial fallback is invoked
    centrally once per accepted step.
  - `thread`: FastBroadcast threading mode used for eligible in-place error-estimate
    broadcasts.

# Keywords

  - `order = Val{1}()`: Taylor expansion order, represented as `Val(N)` for a positive
    integer `N`.
  - `stage_limiter! = trivial_limiter!`: stage-limiter callback with signature
    `limiter!(u, integrator, p, t)`.
  - `step_limiter! = trivial_limiter!`: deprecated algorithm-level fallback with
    signature `limiter!(u, integrator, p, t)`. When no solve-level `step_limiter` is
    supplied, a nontrivial fallback is invoked centrally once per accepted step.
  - `thread = Serial()`: FastBroadcast threading mode. Use a threading mode such as
    `Threaded()` for eligible broadcasts when Julia has multiple threads.

# Returns

An `ExplicitTaylor` algorithm object for use with `solve` or `init`.

# Examples

```julia
using OrdinaryDiffEqTaylorSeries: ExplicitTaylor
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0))
sol = solve(
    prob, ExplicitTaylor(order = Val(8)),
    abstol = 1.0e-10, reltol = 1.0e-10,
)
sol(1.0)
```
"""
Base.@kwdef struct ExplicitTaylor{P, StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    order::Val{P} = Val{1}()
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

"""
    ExplicitTaylorAdaptiveOrder(; min_order = Val{1}(), max_order = Val{10}(),
        stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = Serial())

Explicit Taylor series method that adapts both the expansion order and the step size. At
each accepted step, the controller compares nearby orders and selects the estimated
lowest-work choice for the next step.

The right-hand side must be traceable by Symbolics.jl and must accept TaylorDiff numbers.
The admissible order window is `min_order:(max_order - 1)` because one additional Taylor
coefficient is required to estimate the local error.

# Fields

  - `min_order`: compile-time `Val` containing the lowest admissible step order.
  - `max_order`: compile-time `Val` containing the Taylor jet order. The highest
    admissible step order is one less than this value.
  - `stage_limiter!`: limiter configuration retained by the common explicit-solver
    interface. The current Taylor step does not invoke this callback.
  - `step_limiter!`: deprecated algorithm-level fallback for the accepted-step limiter.
    When no solve-level `step_limiter` is supplied, a nontrivial fallback is invoked
    centrally once per accepted step.
  - `thread`: FastBroadcast threading mode used for eligible in-place error-estimate
    broadcasts.

# Keywords

  - `min_order = Val{1}()`: lowest order considered by the adaptive-order controller,
    represented as `Val(N)`.
  - `max_order = Val{10}()`: Taylor jet order, represented as `Val(N)`. It must be
    greater than `min_order`.
  - `stage_limiter! = trivial_limiter!`: stage-limiter callback with signature
    `limiter!(u, integrator, p, t)`.
  - `step_limiter! = trivial_limiter!`: deprecated algorithm-level fallback with
    signature `limiter!(u, integrator, p, t)`. When no solve-level `step_limiter` is
    supplied, a nontrivial fallback is invoked centrally once per accepted step.
  - `thread = Serial()`: FastBroadcast threading mode. Use a threading mode such as
    `Threaded()` for eligible broadcasts when Julia has multiple threads.

# Returns

An `ExplicitTaylorAdaptiveOrder` algorithm object for use with `solve` or `init`.

# Throws

`solve` throws `ArgumentError` when `max_order <= min_order`.

# Examples

```julia
using OrdinaryDiffEqTaylorSeries: ExplicitTaylorAdaptiveOrder
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 1.0))
sol = solve(
    prob, ExplicitTaylorAdaptiveOrder(min_order = Val(4), max_order = Val(9)),
    abstol = 1.0e-10, reltol = 1.0e-10,
)
sol(1.0)
```
"""
Base.@kwdef struct ExplicitTaylorAdaptiveOrder{P, Q, StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    min_order::Val{P} = Val{1}()
    max_order::Val{Q} = Val{10}()
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

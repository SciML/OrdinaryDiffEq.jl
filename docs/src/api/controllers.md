# Controller API

Starting with OrdinaryDiffEq v7, the adaptive step-size controller is a
first-class object. Each algorithm picks one through
[`default_controller`](@ref), and users override it by passing
`controller = ...` to `solve` / `init`. Tuning knobs like `qmin`, `qmax`,
`gamma`, `qsteady_min/max`, `beta1`, `beta2`, and `failfactor` live on the
controller — not on `integrator.opts` — and the per-step state
(`q11`, `qold`, `dtacc`, `erracc`, …) lives on the controller's
**cache**, which is owned by `integrator.controller_cache`.

This page documents:

1. The built-in controllers users can pick from.
2. The [`CommonControllerOptions`](@ref) bundle of standard step-size
   knobs that every controller composes.
3. The accessor functions the integrator uses to read knobs
   (`get_qmin`, `get_qmax`, …).
4. The minimal contract you implement when adding a new controller.

## Picking a controller

`solve(prob, alg; controller = MyController())` overrides the default
chosen by `default_controller(QT, alg)`. The built-in choices fall into
two families.

### Generic controllers

These contain the full step-size logic themselves and work with any
adaptive algorithm.

```@docs
OrdinaryDiffEqCore.IController
OrdinaryDiffEqCore.PIController
OrdinaryDiffEqCore.PIDController
OrdinaryDiffEqCore.PredictiveController
```

### Algorithm-integrated controllers

For algorithms whose step-size logic is tied to internal state (the
BDF order, the Nordsieck history, extrapolation tableau choice, …),
the controller is a thin wrapper that just owns the knobs; the real
update lives in the algorithm-specific `stepsize_controller!` /
`step_accept_controller!` / `step_reject_controller!` methods.

```@docs
OrdinaryDiffEqBDF.BDFController
OrdinaryDiffEqNordsieck.JVODEController
OrdinaryDiffEqExtrapolation.ExtrapolationController
ImplicitDiscreteSolve.KantorovichTypeController
```

### Composite controllers

```@docs
OrdinaryDiffEqCore.CompositeController
```

## CommonControllerOptions

Every controller embeds a `CommonControllerOptions{T}` as its `basic`
field. This holds the seven standard step-size scalars and is the
target of [`resolve_basic`](@ref).

```@docs
OrdinaryDiffEqCore.CommonControllerOptions
OrdinaryDiffEqCore.resolve_basic
```

## Reading knobs from the integrator

Hot-path code (the integrator's rejection logic, the implicit Newton
loop, the per-step controller logic itself) reads knobs through these
accessor functions rather than directly off the controller struct.
That layer of indirection lets `CompositeControllerCache` delegate to
the currently-active sub-cache.

```@docs
OrdinaryDiffEqCore.get_qmin
OrdinaryDiffEqCore.get_qmax
OrdinaryDiffEqCore.get_qmax_first_step
OrdinaryDiffEqCore.get_gamma
OrdinaryDiffEqCore.get_qsteady_min
OrdinaryDiffEqCore.get_qsteady_max
OrdinaryDiffEqCore.get_failfactor
OrdinaryDiffEqCore.get_current_qmax
OrdinaryDiffEqCore.get_EEst
OrdinaryDiffEqCore.set_EEst!
```

## Per-algorithm defaults

Algorithms specialize these single-argument functions to set their own
default values for the standard step-size knobs. The defaults are
consulted by [`resolve_basic`](@ref) at `setup_controller_cache` time
for any knob the user did not override.

```@docs
OrdinaryDiffEqCore.qmin_default
OrdinaryDiffEqCore.qmax_default
OrdinaryDiffEqCore.qmax_first_step_default
OrdinaryDiffEqCore.gamma_default
OrdinaryDiffEqCore.qsteady_min_default
OrdinaryDiffEqCore.qsteady_max_default
OrdinaryDiffEqCore.failfactor_default
OrdinaryDiffEqCore.beta1_default
OrdinaryDiffEqCore.beta2_default
OrdinaryDiffEqCore.default_controller
```

For example, the BDF family overrides `qmax_default(::QNDF) = 5 // 1`
and `qsteady_max_default(::QNDF) = 12 // 10` to match the historical
`QNDF(qmax=5, qsteady_max=12//10)` defaults.

## Implementing a custom controller

To plug in a new step-size strategy, you implement two types — the
controller (immutable knob holder) and the controller cache (mutable
per-solve state) — and a small set of methods.

### Required types

```julia
struct MyController{B} <: OrdinaryDiffEqCore.AbstractController
    basic::B  # NamedTuple of overrides (unresolved) or CommonControllerOptions{T} (resolved)
end

mutable struct MyControllerCache{T, E} <: OrdinaryDiffEqCore.AbstractControllerCache
    controller::MyController{OrdinaryDiffEqCore.CommonControllerOptions{T}}
    # ... per-step scratch state owned by the controller ...
    EEst::E   # scalar error estimate; required unless you override `get_EEst` / `set_EEst!`
end
```

The cache must expose the scalar error estimate via an `EEst` field, or
provide its own [`get_EEst`](@ref) / [`set_EEst!`](@ref) overrides.

### Required methods

```@docs
OrdinaryDiffEqCore.AbstractController
OrdinaryDiffEqCore.AbstractControllerCache
OrdinaryDiffEqCore.setup_controller_cache
OrdinaryDiffEqCore.stepsize_controller!
OrdinaryDiffEqCore.step_accept_controller!
OrdinaryDiffEqCore.step_reject_controller!
```

In summary:

| Method | What it does | When called |
|---|---|---|
| `setup_controller_cache(alg, alg_cache, controller, ::Type{EEstT})` | Allocate the controller cache, resolve any unresolved [`CommonControllerOptions`](@ref). | `init`, once per solve. |
| `stepsize_controller!(integrator, cache, alg)` | Compute and return the step-size factor `q` based on the current error estimate. | Every step, after `perform_step!`. |
| `step_accept_controller!(integrator, cache, alg, q)` | Return the new `dt` for the next step when the current one is accepted. | Once per accepted step. |
| `step_reject_controller!(integrator, cache, alg)` | Mutate `integrator.dt` so the rejected step can be retried. | Once per rejected step. |

### Optional methods

These have library-provided defaults that work for most controllers but
are commonly overridden for algorithm-integrated controllers (BDF,
JVODE, …) where the algorithm itself owns the step-size logic.

```@docs
OrdinaryDiffEqCore.accept_step_controller
OrdinaryDiffEqCore.post_newton_controller!
OrdinaryDiffEqCore.reinit_controller!
OrdinaryDiffEqCore.sync_controllers!
OrdinaryDiffEqCore.reset_alg_dependent_opts!
```

| Method | Default | When to override |
|---|---|---|
| `accept_step_controller(integrator, cache, alg)::Bool` | `get_EEst(integrator) <= 1` | Controllers that accept on a different criterion (e.g. `PIDController` uses `dt_factor >= accept_safety`). |
| `post_newton_controller!(integrator, cache, alg)` | Shrinks `integrator.dt` by `get_failfactor(integrator)` | Implicit-solver controllers that also reduce the BDF order on Newton failure (`BDFController`, `JVODEController`). |
| `reinit_controller!(integrator, cache)` | no-op | Stateful controllers that need to reset `q11`, `qold`, `dtacc`, etc. on `reinit!`. |
| `sync_controllers!(cache1, cache2)` | no-op | Composite-alg switching, when state must transfer between sub-controllers. |
| `reset_alg_dependent_opts!(cache, alg1, alg2)` | no-op | Re-derive `beta1` / `beta2` etc. when a composite algorithm switches between branches. |
| `get_qmin(integrator, cache)` (and `get_qmax` / `get_gamma` / `get_qsteady_min` / `get_qsteady_max` / `get_qmax_first_step` / `get_failfactor`) | Reads `cache.controller.basic.X` | Controllers that don't embed a [`CommonControllerOptions`](@ref) or store the knobs somewhere other than `cache.controller.basic`. |
| `get_EEst(cache)` / `set_EEst!(cache, val)` | Read/write the `EEst` field on the cache | Caches that store the error differently (vector, per-stage, …). |

### Worked example

Here is a stripped-down implementation of an I-controller using the new
interface. It composes [`CommonControllerOptions`](@ref) so users can
pass `qmin` / `qmax` / `gamma` as kwargs, and it picks up
algorithm-specific defaults automatically.

A note on conventions before the code: the controller protocol carries
`q` as a *divisor* — the new step size is `dt / q`, so `q > 1` shrinks
the step, `q < 1` grows it, and the controller clamps `q` to
`[inv(qmax), inv(qmin)]`. Every built-in controller follows this; new
controllers should too so they compose with each other and with
[`CompositeController`](@ref).

```julia
using OrdinaryDiffEqCore
using OrdinaryDiffEqCore: AbstractController, AbstractControllerCache,
                         CommonControllerOptions, resolve_basic,
                         _resolved_QT, get_EEst, get_current_adaptive_order,
                         get_current_qmax, fastpower
import OrdinaryDiffEqCore: setup_controller_cache, stepsize_controller!,
                          step_accept_controller!, step_reject_controller!

# Type-stability trick: `basic` is either a `NamedTuple` of user-supplied
# overrides (before setup_controller_cache resolves it) or a fully-resolved
# `CommonControllerOptions{T}` (after). The Union constraint lets the
# compiler still specialize on both branches.
struct MyIController{B <: Union{NamedTuple, CommonControllerOptions}} <: AbstractController
    basic::B
end

# Keyword form: the user's `qmin = …`, `qmax = …`, … ride along as a
# `NamedTuple` until setup_controller_cache merges them with the
# algorithm's defaults.
MyIController(; kwargs...) = MyIController(NamedTuple(kwargs))

# The cache holds the resolved controller (concrete `CommonControllerOptions{T}`)
# and any per-solve scratch state — for an I-controller, just the scalar EEst.
mutable struct MyIControllerCache{T, E} <: AbstractControllerCache
    controller::MyIController{CommonControllerOptions{T}}
    EEst::E
end

# Called once at `init` time. Picks the scalar type `QT` for the step-size
# factors, then fills in every unset knob from `qmin_default(alg)`,
# `qmax_default(alg)`, … via `resolve_basic`. After this, every read of
# `cache.controller.basic.X` is type-stable.
function setup_controller_cache(alg, alg_cache, controller::MyIController, ::Type{E}) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT)
    return MyIControllerCache{QT, E}(MyIController(basic), oneunit(E))
end

# Called after each `perform_step!`. Computes the divisor `q` from the
# scaled error estimate using the standard formula
# `q = EEst^(1/(p+1)) / gamma`, clamped to `[inv(qmax), inv(qmin)]`.
# `get_current_qmax` returns a looser `qmax_first_step` on the very first
# step (CVODE-style) and the regular `qmax` afterwards.
@inline function stepsize_controller!(integrator, cache::MyIControllerCache, alg)
    (; qmin, qmax, gamma) = cache.controller.basic
    qmax = get_current_qmax(integrator, qmax)
    EEst = get_EEst(integrator)             # scaled error estimate from perform_step!
    if iszero(EEst)
        return inv(qmax)                    # grow as fast as allowed
    end
    expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
    qtmp = fastpower(EEst, expo) / gamma    # raw factor predicted by I controller
    return max(inv(qmax), min(inv(qmin), qtmp))
end

# Called once per accepted step. Returns the new `dt` for the next step.
# The deadband (`qsteady_min ≤ q ≤ qsteady_max`) holds `dt` constant when
# the proposed change is small — useful for avoiding jitter, and required
# by Jacobian-reusing implicit solvers.
function step_accept_controller!(integrator, cache::MyIControllerCache, alg, q)
    (; qsteady_min, qsteady_max) = cache.controller.basic
    if qsteady_min <= q <= qsteady_max
        q = one(q)                          # within deadband — keep dt
    end
    return integrator.dt / q                # new dt for the next attempt
end

# Called once per rejected step. Must mutate `integrator.dt` so the
# step gets retried with a smaller value. Here we just shrink by `qmax`
# (the most aggressive contraction the controller is allowed to make).
step_reject_controller!(integrator, cache::MyIControllerCache, alg) =
    (integrator.dt = integrator.dt / get_current_qmax(integrator, cache.controller.basic.qmax))
```

You can then pass it through `solve`:

```julia
sol = solve(prob, Tsit5(); controller = MyIController(qmax = 5))
```

`Tsit5()` is non-stiff so it would normally get `PIController`; here we
swap it for the custom one. The `qmax = 5` rides through as part of the
unresolved `NamedTuple`; [`resolve_basic`](@ref) fills in `qmin`,
`gamma`, `qsteady_min`, `qsteady_max`, `qmax_first_step`, and
`failfactor` from `qmin_default(::Tsit5)`, `gamma_default(::Tsit5)`,
etc., at `init` time.

### Adding a controller-specific knob

If your controller relies on a knob that doesn't exist on
[`CommonControllerOptions`](@ref), add it as a new field on your
controller struct (alongside `basic`). The built-in controllers follow
this pattern: `PIController` keeps `beta1` / `beta2` / `qoldinit`
separate from the standard knobs, and `PIDController` keeps `beta` /
`accept_safety` / `limiter`. Separating them keeps
[`CommonControllerOptions`](@ref) small enough to be useful as the
"every adaptive controller has these" surface, while still letting
specialized controllers carry their own tuning state. The new field
moves through `setup_controller_cache` the same way `basic` does —
unresolved on the user-constructed controller, then converted to the
final type `QT` and stored on the resolved controller when the cache
is built.

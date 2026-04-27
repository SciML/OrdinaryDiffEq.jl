# OrdinaryDiffEq.jl v7 Breaking Changes

This release bumps to **SciMLBase v3**, **RecursiveArrayTools v4**, and includes breaking changes across **DiffEqBase**, **OrdinaryDiffEqCore**, and all solver sublibraries.

## Themes of the v7 release

Most of the breaking changes fall into a small set of recurring themes. Keep these in mind while reading the migration table — they explain why an individual change exists and often suggest the right migration direction:

- **Time to first solve (TTFS) reduction.** Direct deps on `Static.jl`, `StaticArrayInterface.jl`, `Polyester.jl`, and `StaticArrays.jl` were dropped; `using OrdinaryDiffEq` now loads only the default solver set; `ODEFunction` switched to `AutoSpecialize`. All of this means less code loaded and more precompilation caching on first `solve`.
- **Type stability everywhere.** All `Bool` solver/`solve` keyword arguments (`autodiff`, `verbose`, `alias`, `lazy`, …) were replaced by typed objects. Passing a `Bool` no longer changes dispatch in ways the compiler cannot specialize on, and the reverse is no longer allowed to silently fall back through slow generic paths.
- **Generality beyond ForwardDiff.** `chunk_size`, `diff_type`, `standardtag`, etc. encoded ForwardDiff-specific or FiniteDiff-specific knobs on every solver. They are replaced by the `ADTypes` interface (`AutoForwardDiff`, `AutoFiniteDiff`, `AutoEnzyme`, `AutoZygote`, …) so every solver automatically generalizes to any AD backend.
- **Controller is now an object, not a pile of `solve` kwargs.** `gamma`, `beta1`, `beta2`, `qmin`, `qmax`, `qsteady_min`, `qsteady_max`, `qoldinit` were moved onto concrete `PIController` / `PIDController` / `IController` / `PredictiveController` structs, and `EEst` moved to the controller cache. This is prep work for pluggable controllers and removes a large amount of dead state from the integrator struct.
- **Cleanup of old re-exports / deprecations.** Functions like `has_destats` (now `has_stats`), `sol.destats` (now `sol.stats`), `DEAlgorithm`/`DEProblem`/`DESolution` abstract types, `tuples()`/`intervals()`, `QuadratureProblem`, `fastpow`, `concrete_solve`, etc. were on a deprecation path for one or more minor releases. v7 removes them.

## Recommended upgrade path

The cleanest path is **not** to jump straight from an old environment onto v7. Most renamed APIs (e.g. `DEAlgorithm` → `AbstractDEAlgorithm`, `u_modified!` → `derivative_discontinuity!`, `has_destats` → `has_stats`, `sol.destats` → `sol.stats`, the `construct*` tableau functions, `alias_u0`/`alias_du0`, `beta1`/`beta2`, PID kwargs) already exist under their new names in **SciMLBase v2** / **OrdinaryDiffEq v6** with deprecation warnings. The recommended sequence is:

1. **Stay on SciMLBase v2 / OrdinaryDiffEq v6.** Update your code to the new names (`has_stats`, `sol.stats`, `AbstractDEAlgorithm`, `derivative_discontinuity!`, `ODEAliasSpecifier`, `ODEVerbosity`, `ADTypes`-based `autodiff`, explicit `controller = …` objects, new tableau names) while the deprecation shims still exist.
2. **Verify your tests pass on v6 with no deprecation warnings.**
3. **Then bump to v7.** At this point your code should compile and run against v7 without further changes aside from the genuinely new breakage (RAT v4 array semantics, ensemble `prob_func`/`output_func` signature, struct type parameter removals, kwargs that truly no longer exist, default changes like `CheckInit` and `williamson_condition=false`).

Doing it in two steps keeps the diff small per step and lets the deprecation warnings on v6 point you at the exact call sites that will break on v7.

## Fallback for RAT v4 indexing

If you cannot update `sol[i]` / `length(sol)` / `eachindex(sol)` call sites yet (see the RAT v4 table below), you can opt back into v3 semantics on a per-solution basis by converting the container type to the ragged variant:

```julia
using RecursiveArrayToolsRaggedArrays
sol_old = RaggedVectorOfArray(sol)   # indexes like v3: sol_old[i] is the i-th timestep
```

`RecursiveArrayToolsRaggedArrays.jl` preserves the previous `AbstractVectorOfArray` indexing behavior (timestep-first, not element-first). This is the escape hatch for code that assumes `sol[i]` returns the i-th timestep. It is, however, recommended that you update to the `sol.u[i]` / `sol[:, i]` style — the ragged wrapper is a compatibility layer, not the canonical API going forward.

---

## RecursiveArrayTools v4

### ODESolution is now an AbstractArray

`AbstractVectorOfArray` (the parent type of `ODESolution`, `RODESolution`, `DAESolution`, etc.) now subtypes `AbstractArray`. This changes the semantics of several common operations:

| Operation | v3 (old) | v4 (new) | Migration |
|---|---|---|---|
| `sol[i]` | Returns i-th timestep (`Vector`) | Returns i-th scalar element (column-major) | Use `sol.u[i]` or `sol[:, i]` |
| `length(sol)` | Number of timesteps | `prod(size(sol))` (total elements) | Use `length(sol.t)` or `length(sol.u)` |
| `eachindex(sol)` | `1:nsteps` | `CartesianIndices(size(sol))` | Use `eachindex(sol.u)` |
| `iterate(sol)` | Iterates over timesteps | Iterates over scalar elements | Use `for u in sol.u` |
| `first(sol)` / `last(sol)` | First/last timestep | First/last scalar element | Use `first(sol.u)` / `last(sol.u)` |
| `map(f, sol)` | Maps over timesteps | Maps over elements | Use `map(f, sol.u)` |
| `maximum(sol)` | Maximum over timesteps | Maximum over all elements | Use `maximum(f, sol.u)` |

**Why:** making `AbstractVectorOfArray <: AbstractArray` lets every generic `AbstractArray` consumer (LinearAlgebra, broadcasting, Zygote adjoints, `StructArrays`, etc.) work on solutions without any special casing in SciMLBase, and deletes a large pile of manual method overrides.

**Migration shortcut:** `sol.u[i]` is the forward-compatible form under both v3 and v4 — if you change every `sol[i]` → `sol.u[i]` in your code base now, it works on both versions.

**Full fallback:** see the "Fallback for RAT v4 indexing" section above — wrapping with `RaggedVectorOfArray` from `RecursiveArrayToolsRaggedArrays.jl` restores v3 indexing.

### EnsembleSolution indexing

The same `AbstractArray` migration applies to `EnsembleSolution` (from `EnsembleProblem` trajectories) and to `EnsembleAnalysis` helpers — an ensemble's `.u` is a `Vector{<:ODESolution}` and the solution itself subtypes `AbstractVectorOfArray`:

| Operation | v3 (old) | v4 (new) | Migration |
|---|---|---|---|
| `sim[j]` | j-th trajectory (`ODESolution`) | j-th scalar element of the flattened container | Use `sim.u[j]` |
| `sim[i, j]` | Trajectory `j`'s i-th timestep (`Matrix` / `Vector`) | Scalar element at column-major position `(i, j)` | Use `sim.u[j].u[i]` |
| `sim[i, j, k]` | Row `k` of trajectory `j`'s i-th timestep | Scalar at `(i, j, k)` | Use `sim.u[j].u[i][k]` |
| `length(sim)` | Number of trajectories | `prod(size(sim))` (total scalar count) | Use `length(sim.u)` |
| `for sol in sim` | Iterate trajectories | Iterate scalar elements (column-major) | Use `for sol in sim.u` |

**Migration shortcut:** as with `ODESolution`, `sim.u[j]` / `sim.u[j].u[i]` / `length(sim.u)` are the forward-compatible forms that work under both v3 and v4.

**EnsembleAnalysis helpers (SciMLBase ≥ 3.4.3):** `get_timestep`, `get_timepoint`, `timeseries_steps_mean/median/quantile/meanvar/meancov/meancor/weighted_meancov`, `EnsembleSummary`, and the weak-error path in `calculate_ensemble_errors` were updated to iterate `sim.u` and use `length(sim.u[1].t)` for the timestep count rather than `length(sim)` / `length(sim.u[1])`. If you pin SciMLBase < 3.4.3 you may hit a `BoundsError` at `calculate_ensemble_errors` on a weak `EnsembleProblem` with `error_estimate`; bump to 3.4.3 or later.

See [SciML/OrdinaryDiffEq.jl#3532](https://github.com/SciML/OrdinaryDiffEq.jl/pull/3532) and [SciML/SciMLBase.jl#1326](https://github.com/SciML/SciMLBase.jl/pull/1326) for the specific v3→v4 ensemble-indexing migrations.

### Other RAT v4 changes

- `zero(VectorOfArray)` now preserves container type (e.g. `StructVector`) via `rewrap`
- Ragged arrays: `size(sol)` reports maximum dimensions; out-of-bounds elements return `zero(T)`
- Removed: custom `length`, `eachindex`, `iterate`, `first`, `last`, `eltype`, `ndims`, `axes`, `any`, `all`, `sum`, `prod`, `mapreduce`, `map` overrides (all inherited from `AbstractArray`)
- Removed: `convert(::Type{AbstractArray}, ...)` (identity since it IS an AbstractArray)
- Removed: `==(::AbstractVectorOfArray, ::AbstractArray)` override

---

## SciMLBase v3

### Renamed APIs

The new names already exist under SciMLBase v2 with deprecation warnings. Update to the v3 names while still on v2 before bumping.

| Old (v2, removed in v3) | New (v2 + v3) | Migration |
|---|---|---|
| `u_modified!(integrator, bool)` | `derivative_discontinuity!(integrator, bool)` | Search-and-replace. `u_modified!` was misleading — the callback system doesn't care whether `u` changed, it cares whether the derivative is discontinuous. |
| `integrator.u_modified` | `integrator.derivative_discontinuity` | Field access rename. |
| `DEAlgorithm` | `AbstractDEAlgorithm` | Abstract type name aligned with SciML-wide `Abstract…` convention. |
| `DEProblem` | `AbstractSciMLProblem` | Same. |
| `DESolution` | `AbstractSciMLSolution` | Same. |
| `sol.destats` | `sol.stats` | Drop the `de` prefix. |
| `has_destats(alg)` | `has_stats(alg)` | Same renaming applied to the trait function. |

### Removed deprecations

All of these printed a deprecation warning under SciMLBase v2. They are gone in v3:

- `has_destats` function → use `has_stats`
- `symbol_to_ReturnCode` and Symbol-to-ReturnCode conversion → use `ReturnCode.*` directly (see next section)
- `syms`/`paramsyms`/`indepsym` kwargs on all `SciMLFunction` constructors → use the problem's `sys` / MTK symbolic interface
- `sol.x` on `AbstractOptimizationSolution` → use `sol.u`
- `prob.lb`/`prob.ub` on `IntegralProblem` → use `prob.domain`
- `sol.minimizer`/`sol.minimum` → use `sol.u` / `sol.objective`
- `tuples()`, `intervals()` iterator functions → iterate over `sol.u` / `sol.t` directly
- `IntegratorTuples`, `IntegratorIntervals`, `TimeChoiceIterator` types
- `QuadratureProblem` alias → use `IntegralProblem`
- `EnsembleProblem` vector-of-problems constructor → use a `prob_func`
- `IntegralProblem` `nout`/`batch` kwargs → set on the integrand function directly
- `SciMLBaseMLStyleExt` extension (`MLStyle` dependency removed)

### `sol.retcode` is a `ReturnCode.T`, not a `Symbol`

Comparing a solution's retcode against a `Symbol` no longer works:

```julia
# Worked on v1/v2, BROKEN on v3
sol.retcode == :Success
```

The `Symbol` return codes were deprecated years ago in favor of the
`ReturnCode.T` enum, with a deprecation warning printed on every use across
the entire v2 series. The deprecation shim has now been removed.

**Migration:** prefer `SciMLBase.successful_retcode(sol)` over equality
against a specific `ReturnCode.*` value. `successful_retcode` correctly
accepts every success-ish return code (`Success`, `StalledSuccess`,
`ExactSolutionLeft`, `ExactSolutionRight`, `FloatingPointLimit`, …), not
just `Success` — so a solver that terminated at an exact solution or hit
a floating-point limit isn't misclassified as a failure.

| Old | New (v3) |
|-----|----------|
| `sol.retcode == :Success` | `SciMLBase.successful_retcode(sol)` |
| `sol.retcode == :Failure` | `!SciMLBase.successful_retcode(sol)` (or match the specific `ReturnCode.Failure` if you really need that exact code) |
| `sol.retcode == :MaxIters` | `sol.retcode == ReturnCode.MaxIters` |
| `sol.retcode == :Default` | `sol.retcode == ReturnCode.Default` |

The full enum is defined in `SciMLBase/src/retcodes.jl` and documented at
<https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#retcodes>.

### Changed defaults

- `ODEFunction{iip}(f)` now uses `AutoSpecialize` by default (was `FullSpecialize`).
  **Why:** drastically more precompilation caching now "works" out of the box — `AutoSpecialize` stores a lightly-specialized method that can be reused across different right-hand side function types instead of baking a fresh specialization per-user. Combined with the default-solver-set reduction (below), this is the single biggest TTFS improvement in v7.
  **Migration:** if you *want* the old behavior (e.g. you are benchmarking and need maximum inlining of `f`), write `ODEFunction{iip, FullSpecialize}(f)` explicitly.
- `is_discrete_time_domain(nothing)` now returns `false` (was `true`). Affects only code that called this trait with `nothing` — replace with an explicit domain.

### Ensemble RNG redesign

- `prob_func(prob, i, repeat)` → `prob_func(prob, ctx)` where `ctx::EnsembleContext`
- `output_func(sol, i)` → `output_func(sol, ctx)`
- New `seed` / `rng` / `rng_func` kwargs on `solve()` for deterministic, thread-count-independent ensemble solves

**Why:** the old `(i, repeat)` signature had no way to plumb a reproducible per-trajectory RNG through, so ensemble results depended on thread count and scheduling. `EnsembleContext` carries an RNG derived from a single top-level seed, so results are reproducible regardless of how many workers run in parallel.

**Migration:**
```julia
# v6
prob_func = (prob, i, repeat) -> remake(prob, u0 = rand(length(prob.u0)))

# v7
prob_func = (prob, ctx) -> remake(prob, u0 = rand(ctx.rng, length(prob.u0)))
```
Within the new signature, `ctx.i` is the trajectory index and `ctx.repeat` is the retry counter if you need them. Pass `solve(ensemble_prob, alg; seed=42, trajectories=N)` for a reproducible run.

### Removed Zygote integer-indexing adjoints for ODESolution

Integer-indexing Zygote adjoints removed; Zygote's generic `AbstractArray` rules handle this under RAT v4. **Migration:** nothing to do if you were using `sol[i]` for AD — the generic path produces the same gradient. Code that explicitly called the `ChainRulesCore`/`Zygote` rule by name should delete those references.

---

## OrdinaryDiffEq v7

### Package scope reduction

`using OrdinaryDiffEq` now loads only the default solver set:

- **Included**: `DefaultODEAlgorithm`, `Tsit5`, `AutoTsit5`, `Vern6`–`Vern9`, `AutoVern6`–`AutoVern9`, `Rosenbrock23`, `Rodas5P`, `FBDF`
- **Not included**: all other solvers require explicit import of their sublibrary

```julia
# v6
using OrdinaryDiffEq
solve(prob, KenCarp4())  # worked

# v7
using OrdinaryDiffEqSDIRK  # explicit import required
solve(prob, KenCarp4())
```

**Why:** TTFS. The old umbrella package loaded every solver family (exponential integrators, symplectic RK, stabilized methods, multirate, Taylor series, …) whether you used them or not. The solver you actually care about, plus its precompilation cache, is now a small fraction of what gets loaded.

**Migration:** add `using OrdinaryDiffEq<Family>` for any non-default solver you use. The sublibrary name is predictable — `KenCarp*`/`TRBDF2` → `OrdinaryDiffEqSDIRK`, `Rosenbrock*`/`Rodas*` except `Rosenbrock23`/`Rodas5P` → `OrdinaryDiffEqRosenbrock`, `RadauIIA*` → `OrdinaryDiffEqFIRK`, etc. Every family has its own `lib/OrdinaryDiffEq<X>` directory in the repo.

### Algorithm struct type parameters removed

All `{CS, AD, FDT, ST, CJ}` type parameters removed from algorithm abstract and concrete types. Over 100 structs affected across BDF, SDIRK, Rosenbrock, FIRK, Extrapolation, ExponentialRK, IMEXMultistep, PDIRK, Newmark, StabilizedIRK, and StochasticDiffEqImplicit.

```julia
# v6
ImplicitEuler{0, AutoForwardDiff{nothing, Nothing}, Nothing, NLNewton{...}, typeof(DEFAULT_PRECS), Val{:forward}(), true, nothing, typeof(trivial_limiter!)}

# v7
ImplicitEuler{AutoForwardDiff{nothing, Nothing}, Nothing, NLNewton{...}, typeof(trivial_limiter!), Nothing}
```

**Why:** the removed parameters (`CS` = chunk size, `FDT` = finite diff type, `ST` = standard tag, `CJ` = concrete Jacobian) are all now carried by the `ADTypes` object on the `autodiff` field, not baked into the algorithm's type. This shrinks the method table, improves compile time, and means a single `autodiff = AutoForwardDiff(chunksize=12)` now communicates what used to take five type parameters.

**Migration:** any code dispatching on `SomeAlg{CS, AD, FDT, ST, CJ}` will break. Most call sites want either `SomeAlg` (no parameters at all) or dispatch on the autodiff type via `alg.autodiff`. If you genuinely need to specialize on chunksize, inspect the `ADTypes` object at runtime.

### Removed keyword arguments from algorithm constructors

Across all implicit/Rosenbrock/BDF/SDIRK/FIRK/Exponential constructors:

| Removed kwarg | Migration |
|---|---|
| `chunk_size` | Set via `autodiff = AutoForwardDiff(chunksize=N)` |
| `diff_type` | Set via `autodiff = AutoFiniteDiff(fdtype=Val(:central))` |
| `standardtag` | Always true; remove the kwarg |
| `precs` | Preconditioners now configured via `linsolve` kwarg (see LinearSolve.jl docs on the `Pl`/`Pr` interface) |
| `controller` (on individual algs) | Set at `solve()` level via `controller = PIController(...)` etc. |

**Why:** `chunk_size` and `diff_type` only meant anything to ForwardDiff and FiniteDiff respectively. Hoisting them onto the `ADTypes` object means every solver now works with **any** AD backend (Enzyme, Zygote, ReverseDiff, Mooncake, …) by just swapping `autodiff=AutoEnzyme()` — no per-solver kwarg surface needed. `standardtag` was always the right default. `precs` moved under `linsolve` because LinearSolve now owns the preconditioner abstraction.

### autodiff: Bool no longer accepted

```julia
# v6
Rosenbrock23(autodiff=true)   # worked
Rosenbrock23(autodiff=false)  # worked

# v7 — must use ADTypes
Rosenbrock23(autodiff=AutoForwardDiff())
Rosenbrock23(autodiff=AutoFiniteDiff())
```

**Why:** type stability. `autodiff::Bool` forced a runtime branch between AD backends that the compiler couldn't specialize through; `autodiff::AbstractADType` dispatches at compile time. This is also what enables the "generalizes beyond ForwardDiff" story — you can now pass `AutoEnzyme()`, `AutoZygote()`, `AutoMooncake()`, etc., and the solver specializes correctly.

### verbose: Bool no longer accepted

```julia
# v6
solve(prob, alg, verbose=false)

# v7 — must use ODEVerbosity
solve(prob, alg, verbose=ODEVerbosity(SciMLLogging.None()))
```

**Why:** same type-stability reason, plus `ODEVerbosity` exposes fine-grained control (separate levels for nonlinear solver, linear solver, initialization, etc.) that a single `Bool` can't express.

### alias: Bool no longer accepted

```julia
# v6
solve(prob, alg, alias=true)

# v7 — must use ODEAliasSpecifier
solve(prob, alg, alias=ODEAliasSpecifier(alias_u0=true))
```

Deprecated `alias_u0` / `alias_du0` keyword shortcuts also removed — they already printed warnings on v6. Update to `alias = ODEAliasSpecifier(alias_u0=…, alias_du0=…)` on v6 first, then upgrade.

**Why:** type stability + finer control. `alias=true` used to mean "alias everything that can be aliased," which was ambiguous when some things could be safely aliased and others couldn't. `ODEAliasSpecifier` makes each aliasable buffer an explicit opt-in.

### Default DAE initialization changed to CheckInit

```julia
# v6 — inconsistent initial conditions silently fixed
solve(prob, Rodas5P())

# v7 — errors on inconsistent initial conditions
solve(prob, Rodas5P())  # throws if u0 is inconsistent

# v7 — explicit opt-in to automatic fixing
solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())
```

**Why:** silently fixing an inconsistent DAE initial condition produced wrong answers when the user's `u0` was actually correct but a modeling bug elsewhere made the system look inconsistent. The new default errors loudly; users who want the old "just fix it" behavior opt in explicitly via `initializealg=BrownFullBasicInit()`.

### Controller refactor

PID controller parameters removed from `solve()`/`init()`:

| Removed kwarg | Migration |
|---|---|
| `gamma` | Pass `controller = PIController(gamma=…)` or `PIDController(...)` |
| `beta1`, `beta2` | Constructor args on `PIController` / `PIDController` |
| `qmin`, `qmax` | Same — on the controller object |
| `qsteady_min`, `qsteady_max` | Same |
| `qoldinit` | Same |

**EEst moved from integrator to controller cache:**

```julia
# v6
integrator.EEst

# v7
OrdinaryDiffEqCore.get_EEst(integrator)
OrdinaryDiffEqCore.set_EEst!(integrator, val)
```

Fields `EEst`, `qold`, `q11`, `erracc`, `dtacc` removed from `ODEIntegrator` struct.

**Why:** controllers are now real, pluggable objects rather than a collection of loose numeric knobs on `solve`. This makes it possible to write a custom controller and just pass `controller = MyController(…)` without adding yet more kwargs to `solve`, and keeps controller state (including `EEst`) inside the controller where it belongs instead of on the integrator.

**Migration:** if you only ever set `gamma`, just write `controller = PIController(gamma=0.9)` (or whatever value). The default controller per algorithm is unchanged; the kwargs are what moved.

### lazy keyword: Bool no longer accepted

`BS5`, `Vern6`–`Vern9`: `lazy` must be `Val{true}()` or `Val{false}()`, not `true`/`false`. **Why:** type stability — same story as `autodiff` / `verbose` / `alias`.

### williamson_condition default changed

All 2N low-storage RK methods: default changed from `williamson_condition=true` to `williamson_condition=false`. **Why:** this optimization only works for mutable `Array`-style state. Having it on by default silently made the method wrong (or errored) for `StaticArrays`, GPU arrays, `ComponentArrays`, etc. Off-by-default is the safe choice; opt in with `williamson_condition=true` when you know your state is a plain `Array`.

### Threading interface changed

| Old (v6) | New (v7) |
|---|---|
| `OrdinaryDiffEq.False()` / `Static.False()` | `Serial()` (from FastBroadcast) |
| `OrdinaryDiffEq.True()` / `Static.True()` | `Threaded()` (from FastBroadcast) |

All `calculate_residuals!` thread argument types changed from `Union{False, True}` to `Union{Serial, Threaded}`.

**Why:** this change is what makes the `Static.jl` dependency removal possible, and `Static.jl` removal is one of the largest single contributions to v7's TTFS reduction. The `Serial`/`Threaded` types from `FastBroadcast` are semantically identical for this purpose and FastBroadcast is already a transitive dep. **Migration:** if your code pattern-matched on `Static.True`/`Static.False` for threading dispatch, swap to `FastBroadcast.Threaded`/`FastBroadcast.Serial`.

---

## Removed dependencies

All four removals are driven by TTFS.

| Package | Replacement | Why |
|---|---|---|
| **Static.jl** | `FastBroadcast.Serial` / `FastBroadcast.Threaded` | Removing Static slashes load time and eliminates a wide surface of `StaticInt`/`StaticBool` specialization that the compiler was re-running for every user |
| **StaticArrayInterface.jl** | `ArrayInterface.ismutable` | The only thing OrdinaryDiffEq actually used from StaticArrayInterface was the mutability query, which ArrayInterface already provides |
| **Polyester.jl** (direct dep) | Moved to weak dep `OrdinaryDiffEqCorePolyesterExt`; requires `using Polyester` to activate | Polyester loads a nontrivial amount of threading infrastructure. Users who don't enable Polyester-threaded solvers no longer pay for it |
| **StaticArrays.jl** (direct dep) | `SVector`/`MVector` in tableaus replaced with `NTuple`/`Vector`; SA moved to test-only | Loading StaticArrays forces compilation of a large generated-function surface for every solver that mentions an `SVector`. Tableaus used them for constants, which `NTuple` expresses just as statically |

---

## DiffEqBase changes

DiffEqBase is now a sublibrary under `lib/DiffEqBase` (migrated from standalone package). Key changes:

- `has_destats` re-export removed (function removed from SciMLBase v3 — use `has_stats`)
- `RECOMPILE_BY_DEFAULT` re-export removed (unused)
- `fastpow` deprecation removed — use `FastPower.fastpower`
- `DEStats` deprecation removed — use `SciMLBase.DEStats`
- `concrete_solve` deprecation removed — use `solve` directly
- `FunctionWrappersWrapper` wrapping updated for new LinearSolve precs interface

**Why migrate DiffEqBase into the monorepo:** it's tightly coupled to OrdinaryDiffEq's internals and most features were only available through OrdinaryDiffEq; keeping them in lockstep releases eliminates a compatibility-bound class of bugs.

---

## Callback changes

### `VectorContinuousCallback` now fires every simultaneous event

Previously, when several conditions of a single `VectorContinuousCallback` crossed zero on the same step, only the first crossing's affect was applied. Other simultaneous events were silently dropped, even though their root was within the step's resolution. Bouncing-balls, multi-contact mechanics, and any callback used as a friction/threshold state machine were all affected.

In v7, `VectorContinuousCallback` resolves all simultaneous events on the same step and dispatches them in one call to the user's `affect!`. See [SciML/OrdinaryDiffEq.jl#3230](https://github.com/SciML/OrdinaryDiffEq.jl/pull/3230) (rolled into the v7 merge [#3242](https://github.com/SciML/OrdinaryDiffEq.jl/pull/3242)) and the follow-up fix [#3549](https://github.com/SciML/OrdinaryDiffEq.jl/pull/3549).

### Breaking: `VectorContinuousCallback` `affect!` signature changed

The `affect!` callback for `VectorContinuousCallback` used to be invoked once per triggering condition, with the index of that single event:

```julia
# v6
affect!(integrator, event_index::Int)
```

In v7 it is invoked **once per step** with a `Vector{Int8}` mask describing every condition's status:

```julia
# v7
affect!(integrator, simultaneous_events::Vector{Int8})
```

Each entry of `simultaneous_events` encodes both whether condition `i` triggered and which way it crossed:

| value | meaning |
|---|---|
| `0`  | condition did not trigger this step |
| `-1` | upcrossing (condition went from negative to positive) |
| `+1` | downcrossing (condition went from positive to negative) |

The vector's length is the callback's `len`; the entries are stable across steps. `affect_neg!` is no longer called for `VectorContinuousCallback` — your single `affect!` handles both crossing directions by inspecting the sign of each nonzero entry.

**Migration:**

```julia
# v6: branch on the event index, optionally a separate affect_neg!
function affect!(integrator, event_index)
    if event_index == 1
        # ball 1 hit the ground
    elseif event_index == 2
        # ball 2 hit the ground
    end
end
cb = VectorContinuousCallback(condition, affect!, affect_neg!, 2)

# v7: branch on which entries fired, sign tells you the direction
function affect!(integrator, simultaneous_events)
    for i in eachindex(simultaneous_events)
        s = simultaneous_events[i]
        s == 0 && continue
        if i == 1
            # ball 1 hit the ground; s == -1 upcrossing, s == +1 downcrossing
        elseif i == 2
            # ball 2 hit the ground
        end
    end
end
cb = VectorContinuousCallback(condition, affect!, 2)
```

The previous "call `affect!` once per event index" behavior was load-bearing for very few users, since the pre-v7 implementation already only invoked it for the *first* simultaneous event anyway. The new shape makes the simultaneous case representable without additional API surface.

**Why a `Vector{Int8}` rather than `Vector{Bool}`:** the sign carries the crossing direction, which previously required a separate `affect_neg!` callback. Folding both into the mask removes the `affect!` / `affect_neg!` split for `VectorContinuousCallback`, so users no longer have to maintain two parallel functions to get up- and down-crossing handling. (The `ContinuousCallback` `affect!` / `affect_neg!` split is unchanged.)

---

## OrdinaryDiffEqCore changes

- `DEOptions` struct: removed `gamma`, `qmax`, `qmin`, `qsteady_max`, `qsteady_min`, `qoldinit`, `controller` fields and the `Controller` type parameter (moved to controller object; see Controller refactor above)
- `ODEIntegrator` struct: removed `EEst`, `qold`, `q11`, `erracc`, `dtacc`, `q` fields and `EEstT`, `QT` type parameters; added `controller_cache` field. **Migration:** `integrator.EEst` → `OrdinaryDiffEqCore.get_EEst(integrator)`; the rest live on `integrator.controller_cache`
- `u_modified` field → `derivative_discontinuity`
- `DEFAULT_PRECS` removed (preconditioners now configured via `linsolve`)
- `ispredictive` / `isstandard` controller traits removed — algorithms now override `default_controller` directly. **Migration:** if you defined a custom algorithm that set `isstandard(::MyAlg) = true`, instead define `default_controller(::MyAlg, args...) = IController(…)`
- `has_chunksize` removed (dead code — chunksize is on the ADTypes object now)
- Backwards-compat positional `Alg(stage_limiter!, step_limiter!)` constructors removed from every explicit RK sublibrary — 99 constructors total across LowStorageRK (44), SSPRK (18), LowOrderRK (22), Verner (5: Vern6/7/8/9, RKV76IIa), HighOrderRK (4: TanYam7, TsitPap8, DP8, PFRK87), SIMDRK (3: MER5v2, MER6v2, RK6v4), QPRK (QPRK98), TaylorSeries (ExplicitTaylor2), and Tsit5. **Migration:** use the keyword form `Alg(; stage_limiter! = my_limiter, step_limiter! = my_limiter)`. The kwarg constructors (and no-arg `Alg()` defaulting both limiters to `trivial_limiter!`) are unchanged.

---

## Tableau changes

Tableau functions in `OrdinaryDiffEqExplicitTableaus` / `OrdinaryDiffEqImplicitTableaus` dropped the `construct` prefix:

```julia
# v6
constructDormandPrince()

# v7
OrdinaryDiffEqExplicitTableaus.DormandPrince()
```

~80+ functions renamed. Functions are no longer exported — qualify explicitly with the sublibrary name. The new names also exist on v6 as deprecations, so you can rename first and bump later.

### DiffEqDevTools tableau constructors removed (DiffEqDevTools v2 → v3, breaking)

Previously, DiffEqDevTools re-exported 105 `construct*` tableau constructors (`constructEuler`, `constructKutta3`, `constructRK4`, `constructRK438Rule`, `constructSSPRK22`, `constructImplicitEuler`, `constructMidpointRule`, `constructTrapezoidalRule`, `constructLobattoIIIA4`, `constructGL2`, …) defined in its own `src/ode_tableaus.jl`. These are **removed** in DiffEqDevTools v3 — the authoritative tableau definitions now live exclusively in `OrdinaryDiffEqExplicitTableaus` / `OrdinaryDiffEqImplicitTableaus` under their renamed bare forms.

**Migration:**

```julia
# v6 (DiffEqDevTools v2)
using DiffEqDevTools
tab = constructRK4()

# v7 (DiffEqDevTools v3)
using OrdinaryDiffEqExplicitTableaus
tab = OrdinaryDiffEqExplicitTableaus.RK4()
```

`DiffEqDevTools.deduce_Butcher_tableau(alg)` (which recovers A/b/c from a live solver) is kept unchanged, as is the `ODERKTableau` / `ExplicitRKTableau` / `ImplicitRKTableau` type surface.

---

## StochasticDiffEq / DelayDiffEq

Same theme as the ODE side — all of the following break identically to the ODE versions, with identical migrations:

- Same `verbose`, `alias`, `autodiff` Bool-to-typed-object changes as ODE solvers
- `alias_u0` / `alias_jumps` / `alias_noise` deprecated kwargs removed → use `ODEAliasSpecifier` (for SDE: `SDEAliasSpecifier`)
- `beta1` / `beta2` PID parameter deprecation warnings removed → set on the controller object
- `initial_order` deprecated kwarg removed from DelayDiffEq
- `ispredictive` / `isstandard` trait definitions removed; algorithms now override `default_controller` directly
- `@static if isdefined` version guards removed (always true on v7)

### StochasticDelayDiffEq deprecated

`StochasticDelayDiffEq.jl` is deprecated. Use `DelayDiffEq.jl` directly — it has supported SDDE problems for some time, and the separate `StochasticDelayDiffEq` wrapper is no longer being maintained. It will not receive a v7-compatible release.

**Migration:** replace `using StochasticDelayDiffEq` with `using DelayDiffEq` (plus `using StochasticDiffEq` if you were relying on the SDE algorithm re-exports). `MethodOfSteps(alg)` and the `SDDEProblem` constructor continue to work from `DelayDiffEq` / `SciMLBase` respectively.

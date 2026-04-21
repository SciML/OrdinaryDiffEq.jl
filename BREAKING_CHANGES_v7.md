# OrdinaryDiffEq.jl v7 Breaking Changes

This release bumps to **SciMLBase v3**, **RecursiveArrayTools v4**, and includes breaking changes across **DiffEqBase**, **OrdinaryDiffEqCore**, and all solver sublibraries.

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

### Other RAT v4 changes

- `zero(VectorOfArray)` now preserves container type (e.g. `StructVector`) via `rewrap`
- Ragged arrays: `size(sol)` reports maximum dimensions; out-of-bounds elements return `zero(T)`
- Removed: custom `length`, `eachindex`, `iterate`, `first`, `last`, `eltype`, `ndims`, `axes`, `any`, `all`, `sum`, `prod`, `mapreduce`, `map` overrides (all inherited from `AbstractArray`)
- Removed: `convert(::Type{AbstractArray}, ...)` (identity since it IS an AbstractArray)
- Removed: `==(::AbstractVectorOfArray, ::AbstractArray)` override

---

## SciMLBase v3

### Renamed APIs

| Old (v2) | New (v3) |
|---|---|
| `u_modified!(integrator, bool)` | `derivative_discontinuity!(integrator, bool)` |
| `integrator.u_modified` | `integrator.derivative_discontinuity` |
| `DEAlgorithm` | `AbstractDEAlgorithm` |
| `DEProblem` | `AbstractSciMLProblem` |
| `DESolution` | `AbstractSciMLSolution` |
| `sol.destats` | `sol.stats` |

### Removed deprecations

- `has_destats` function
- `symbol_to_ReturnCode` and Symbol-to-ReturnCode conversion
- `syms`/`paramsyms`/`indepsym` kwargs from all `SciMLFunction` constructors
- `sol.x` on `AbstractOptimizationSolution` (use `sol.u`)
- `prob.lb`/`prob.ub` on `IntegralProblem` (use `prob.domain`)
- `sol.minimizer`/`sol.minimum`
- `tuples()`, `intervals()` iterator functions
- `IntegratorTuples`, `IntegratorIntervals`, `TimeChoiceIterator` types
- `QuadratureProblem` alias
- `EnsembleProblem` vector-of-problems constructor
- `IntegralProblem` `nout`/`batch` kwargs
- `SciMLBaseMLStyleExt` extension (MLStyle dependency removed)

### Changed defaults

- `ODEFunction{iip}(f)` now uses `AutoSpecialize` by default (was `FullSpecialize`)
- `is_discrete_time_domain(nothing)` now returns `false` (was `true`)

### Ensemble RNG redesign

- `prob_func(prob, i, repeat)` → `prob_func(prob, ctx)` where `ctx::EnsembleContext`
- `output_func(sol, i)` → `output_func(sol, ctx)`
- New `seed`/`rng`/`rng_func` kwargs on `solve()` for deterministic, thread-count-independent ensemble solves

### Removed Zygote integer-indexing adjoints for ODESolution

Integer-indexing Zygote adjoints removed; Zygote's generic `AbstractArray` rules handle this under RAT v4.

---

## OrdinaryDiffEq v7

### Package scope reduction

`using OrdinaryDiffEq` now loads only the default solver set:

- **Included**: `DefaultODEAlgorithm`, `Tsit5`, `AutoTsit5`, `Vern6`–`Vern9`, `AutoVern6`–`AutoVern9`, `Rosenbrock23`, `Rodas5P`, `FBDF`
- **Not included**: All other solvers require explicit import of their sublibrary

```julia
# v6
using OrdinaryDiffEq
solve(prob, KenCarp4())  # worked

# v7
using OrdinaryDiffEqSDIRK  # explicit import required
solve(prob, KenCarp4())
```

### Algorithm struct type parameters removed

All `{CS, AD, FDT, ST, CJ}` type parameters removed from algorithm abstract and concrete types. Over 100 structs affected across BDF, SDIRK, Rosenbrock, FIRK, Extrapolation, ExponentialRK, IMEXMultistep, PDIRK, Newmark, StabilizedIRK, and StochasticDiffEqImplicit.

```julia
# v6
ImplicitEuler{0, AutoForwardDiff{nothing, Nothing}, Nothing, NLNewton{...}, typeof(DEFAULT_PRECS), Val{:forward}(), true, nothing, typeof(trivial_limiter!)}

# v7
ImplicitEuler{AutoForwardDiff{nothing, Nothing}, Nothing, NLNewton{...}, typeof(trivial_limiter!), Nothing}
```

Any code dispatching on `SomeAlg{CS, AD, FDT, ST, CJ}` will break.

### Removed keyword arguments from algorithm constructors

Across all implicit/Rosenbrock/BDF/SDIRK/FIRK/Exponential constructors:

| Removed kwarg | Migration |
|---|---|
| `chunk_size` | Set via `autodiff = AutoForwardDiff(chunksize=N)` |
| `diff_type` | Set via `autodiff = AutoFiniteDiff(fdtype=Val(:central))` |
| `standardtag` | Always true; removed |
| `precs` | Preconditioners now configured via `linsolve` kwarg |
| `controller` (on individual algs) | Set at `solve()` level via `controller` kwarg |

### autodiff: Bool no longer accepted

```julia
# v6
Rosenbrock23(autodiff=true)   # worked
Rosenbrock23(autodiff=false)  # worked

# v7 — must use ADTypes
Rosenbrock23(autodiff=AutoForwardDiff())
Rosenbrock23(autodiff=AutoFiniteDiff())
```

### verbose: Bool no longer accepted

```julia
# v6
solve(prob, alg, verbose=false)

# v7 — must use ODEVerbosity
solve(prob, alg, verbose=ODEVerbosity(SciMLLogging.None()))
```

### alias: Bool no longer accepted

```julia
# v6
solve(prob, alg, alias=true)

# v7 — must use ODEAliasSpecifier
solve(prob, alg, alias=ODEAliasSpecifier(alias_u0=true))
```

Deprecated `alias_u0`/`alias_du0` keyword shortcuts also removed.

### Default DAE initialization changed to CheckInit

```julia
# v6 — inconsistent initial conditions silently fixed
solve(prob, Rodas5P())

# v7 — errors on inconsistent initial conditions
solve(prob, Rodas5P())  # throws if u0 is inconsistent

# v7 — explicit opt-in to automatic fixing
solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())
```

### Controller refactor

PID controller parameters removed from `solve()`/`init()`:

| Removed kwarg | Migration |
|---|---|
| `gamma` | Pass `controller = PIController(...)` or `PIDController(...)` |
| `beta1`, `beta2` | Part of controller object |
| `qmin`, `qmax` | Part of controller object |
| `qsteady_min`, `qsteady_max` | Part of controller object |
| `qoldinit` | Part of controller object |

**EEst moved from integrator to controller cache:**

```julia
# v6
integrator.EEst

# v7
OrdinaryDiffEqCore.get_EEst(integrator)
OrdinaryDiffEqCore.set_EEst!(integrator, val)
```

Fields `EEst`, `qold`, `q11`, `erracc`, `dtacc` removed from `ODEIntegrator` struct.

### lazy keyword: Bool no longer accepted

`BS5`, `Vern6`–`Vern9`: `lazy` must be `Val{true}()` or `Val{false}()`, not `true`/`false`.

### williamson_condition default changed

All 2N low-storage RK methods: default changed from `williamson_condition=true` to `williamson_condition=false`. This optimization only works for `Array` types and is now opt-in.

### Threading interface changed

| Old (v6) | New (v7) |
|---|---|
| `OrdinaryDiffEq.False()` / `Static.False()` | `Serial()` (from FastBroadcast) |
| `OrdinaryDiffEq.True()` / `Static.True()` | `Threaded()` (from FastBroadcast) |

All `calculate_residuals!` thread argument types changed from `Union{False, True}` to `Union{Serial, Threaded}`.

---

## Removed dependencies

| Package | Replacement |
|---|---|
| **Static.jl** | `FastBroadcast.Serial` / `FastBroadcast.Threaded` |
| **StaticArrayInterface.jl** | `ArrayInterface.ismutable` |
| **Polyester.jl** (direct dep) | Moved to weak dep `OrdinaryDiffEqCorePolyesterExt`; requires `using Polyester` to activate |
| **StaticArrays.jl** (direct dep) | `SVector`/`MVector` in tableaus replaced with `NTuple`/`Vector`; SA moved to test-only |

---

## DiffEqBase changes

DiffEqBase is now a sublibrary under `lib/DiffEqBase` (migrated from standalone package). Key changes:

- `has_destats` re-export removed (function removed from SciMLBase v3)
- `RECOMPILE_BY_DEFAULT` re-export removed (unused)
- `fastpow` deprecation removed (use `FastPower.fastpower`)
- `DEStats` deprecation removed (use `SciMLBase.DEStats`)
- `concrete_solve` deprecation removed
- `FunctionWrappersWrapper` wrapping updated for new LinearSolve precs interface

---

## OrdinaryDiffEqCore changes

- `DEOptions` struct: removed `gamma`, `qmax`, `qmin`, `qsteady_max`, `qsteady_min`, `qoldinit`, `controller` fields and the `Controller` type parameter
- `ODEIntegrator` struct: removed `EEst`, `qold`, `q11`, `erracc`, `dtacc`, `q` fields and `EEstT`, `QT` type parameters; added `controller_cache` field
- `u_modified` field → `derivative_discontinuity`
- `DEFAULT_PRECS` removed
- `ispredictive`/`isstandard` controller traits removed — algorithms now override `default_controller` directly
- `has_chunksize` removed (dead code)
- Backwards-compat positional constructors removed from LowStorageRK (44 constructors) and Tsit5

---

## Tableau changes

Tableau functions in `OrdinaryDiffEqExplicitTableaus` / `OrdinaryDiffEqImplicitTableaus` dropped the `construct` prefix:

```julia
# v6
constructDormandPrince()

# v7
OrdinaryDiffEqExplicitTableaus.DormandPrince()
```

~80+ functions renamed. Functions are no longer exported.

---

## StochasticDiffEq / DelayDiffEq

- Same `verbose`, `alias`, `autodiff` changes as ODE solvers
- `alias_u0`/`alias_jumps`/`alias_noise` deprecated kwargs removed
- `beta1`/`beta2` PID parameter deprecation warnings removed
- `initial_order` deprecated kwarg removed from DelayDiffEq
- SDDE (StochasticDDE) code paths removed from DelayDiffEq
- `ispredictive`/`isstandard` trait definitions removed; direct `default_controller` overrides used
- `@static if isdefined` version guards removed (always true on v7)

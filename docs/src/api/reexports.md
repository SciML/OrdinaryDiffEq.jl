# Re-exported SciMLBase API

Every OrdinaryDiffEq solver sublibrary (`OrdinaryDiffEqTsit5`, `OrdinaryDiffEqFeagin`, …)
re-exports part of [SciMLBase](https://github.com/SciML/SciMLBase.jl), so that
`using OrdinaryDiffEqFeagin` alone is enough to write and run a solve, exactly as
`using DifferentialEquations` is. This page defines *which* names that is, and why.

Sublibraries that name their re-exports use exactly the list below. The remaining
sublibraries still re-export all of SciMLBase, which is a superset of this list, and are
migrated to it separately so that no release both narrows exports and changes behaviour.

## The rule

A SciMLBase name is re-exported if and only if it belongs to the **user-facing common
interface** — the part of SciMLBase a user writes in a script. Concretely, a name is
re-exported when it falls in one of these five closed categories:

1. **Problem and function types** — the concrete SciMLBase types in the
   [Common Interface API](@ref). This is the set OrdinaryDiffEq uses for ODEs, split ODEs,
   second-order and dynamical ODEs, DAEs, discrete maps, and ensembles. Problem domains solved
   by other libraries, such as optimization, nonlinear equations, SDEs, DDEs, and boundary
   value problems, are not re-exported.
2. **Solving and solution inspection** — `solve`, `solve!`, `init`, `step!`, `remake`,
   `ReturnCode`, `successful_retcode`, the statistics types, and the specialization and
   initialization options accepted by problem constructors and `solve`.
3. **Callbacks** — the callback types and their root-finding options.
4. **Integrator interface** — the mutating and querying verbs documented for use *inside* a
   callback or an `init`/`step!` loop (`u_modified!`, `add_tstop!`, `set_u!`, `get_du`, …).
5. **Ensembles** — the ensemble algorithms and `EnsembleAnalysis`.

Everything else SciMLBase exports stays behind the `SciMLBase.` qualifier. Those exclusions
are categories, not case-by-case decisions:

  - **Abstract types** (`AbstractODEProblem`, `AbstractDEAlgorithm`, …) — dispatch hooks for
    package authors, not something a user constructs.
  - **Solver-author API** — trait predicates (`isinplace`, `has_jac`, `allowscomplex`,
    `requiresgradient`, …), construction helpers (`build_solution`), wrapper and cache types,
    and macros. These are the interface a *new solver package* implements, not the interface a
    user calls.
  - **Internal dispatch tags** (`StandardODEProblem`, `…AliasSpecifier` internals, originator
    types).
  - **Other problem domains** (`OptimizationProblem`, `NonlinearProblem`, `SDEProblem`, …) —
    their solver packages own those interfaces.

Two consequences worth stating explicitly, because they are the point of having a rule:

  - **The list follows real usage.** It is not "whatever this package happened to export
    before"; the problem and function types come from OrdinaryDiffEq's documented common
    interface.
  - **The list is uniform.** Every solver sublibrary re-exports the same names. A name is never
    available from one sublibrary and missing from another, so moving a script between solver
    packages never produces an `UndefVarError`.

The rule is enforced, not just described: `test/qa/reexport_tests.jl` reads the Common
Interface API and asserts that its concrete SciMLBase problem and function types exactly match
this page. It also checks that the enumerated names are exported by SciMLBase and that every
selective sublibrary uses exactly the same list.

## Do not rely on `using` for solver names

Only the names below come from SciMLBase. Algorithm names (`Tsit5`, `Feagin14`, `FBDF`, …) come
from the solver sublibrary that defines them, and everything else in SciMLBase remains reachable
as `SciMLBase.name`.

## The re-exported names

### Problem and function types

```julia
DAEProblem, DiscreteProblem, DynamicalODEProblem, EnsembleProblem, ODEProblem,
SecondOrderODEProblem, SplitODEProblem, DAEFunction, DiscreteFunction,
DynamicalODEFunction, ODEFunction, SplitFunction
```

### Solving and solution inspection

```julia
solve, solve!, init, step!,
remake, ReturnCode, successful_retcode, DEStats,
NLStats, NullParameters, AutoSpecialize, FullSpecialize,
NoSpecialize, FunctionWrapperSpecialize, CheckInit, NoInit,
OverrideInit
```

### Callbacks

```julia
CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
LeftRootFind, RightRootFind, NoRootFind
```

### Integrator interface

```julia
add_saveat!, add_tstop!, auto_dt_reset!, change_t_via_interpolation!,
check_error, check_error!, first_tstop, get_dt,
get_du, get_du!, get_proposed_dt, get_tmp_cache,
pop_tstop!, reinit!, savevalues!, set_abstol!,
set_proposed_dt!, set_reltol!, set_t!, set_u!,
set_ut!, terminate!, u_modified!
```

### Ensembles

```julia
EnsembleAnalysis, EnsembleDistributed, EnsembleSerial, EnsembleSplitThreads,
EnsembleSummary, EnsembleThreads
```

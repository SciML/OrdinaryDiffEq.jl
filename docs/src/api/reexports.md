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

1. **Problem and function types** — every concrete `…Problem` and `…Function` type SciMLBase
   exports. Membership is mechanical: the name ends in `Problem`/`Function`, the type is not
   abstract, and it is not a dispatch tag (`Abstract…`, `Standard…`). Every problem that can be
   *stated* through DifferentialEquations.jl can therefore be stated through any solver
   sublibrary, whether or not that particular sublibrary can solve it.
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

Two consequences worth stating explicitly, because they are the point of having a rule:

  - **The list does not depend on history.** It is not "whatever this package happened to
    export before"; a name is in or out purely by category. Re-deriving the list from a fresh
    SciMLBase gives the same answer.
  - **The list is uniform.** Every solver sublibrary re-exports the same names. A name is never
    available from one sublibrary and missing from another, so moving a script between solver
    packages never produces an `UndefVarError`.

The rule is enforced, not just described: `test/qa/qa_tests.jl` re-derives categories 1 from a
live SciMLBase, checks the enumerated names in categories 2–5 are still exported by SciMLBase,
and asserts that all sublibraries and this page list exactly the same names. Adding a problem
type to SciMLBase therefore fails CI until the sublibraries re-export it.

## Do not rely on `using` for solver names

Only the names below come from SciMLBase. Algorithm names (`Tsit5`, `Feagin14`, `FBDF`, …) come
from the solver sublibrary that defines them, and everything else in SciMLBase remains reachable
as `SciMLBase.name`.

## The re-exported names

### Problem and function types

```julia
AnalyticalProblem, BVProblem, ConvexOptimizationProblem, DAEProblem,
DDEProblem, DiscreteProblem, DynamicalDDEProblem, DynamicalODEProblem,
DynamicalSDEProblem, EigenvalueProblem, EnsembleProblem, HomotopyProblem,
ImmutableNonlinearProblem, ImmutableODEProblem, ImplicitDiscreteProblem, IncrementingODEProblem,
IntegralProblem, IntervalNonlinearProblem, LinearProblem, NoiseProblem,
NonlinearLeastSquaresProblem, NonlinearProblem, ODEProblem, OptimizationProblem,
PDEProblem, RODEProblem, SCCNonlinearProblem, SDDEProblem,
SDEProblem, SampledIntegralProblem, SecondOrderBVProblem, SecondOrderDDEProblem,
SecondOrderODEProblem, SplitODEProblem, SplitSDEProblem, SteadyStateProblem,
TwoPointBVProblem, TwoPointSecondOrderBVProblem, WeightedEnsembleProblem, BVPFunction,
BatchIntegralFunction, DAEFunction, DDEFunction, DiscreteFunction,
DynamicalBVPFunction, DynamicalDDEFunction, DynamicalODEFunction, DynamicalSDEFunction,
HomotopyNonlinearFunction, ImplicitDiscreteFunction, IncrementingODEFunction, IntegralFunction,
IntervalNonlinearFunction, MultiObjectiveOptimizationFunction, NonlinearFunction, ODEFunction,
ODEInputFunction, OptimizationFunction, RODEFunction, SDDEFunction,
SDEFunction, SplitFunction, SplitSDEFunction, TwoPointBVPFunction,
TwoPointDynamicalBVPFunction
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


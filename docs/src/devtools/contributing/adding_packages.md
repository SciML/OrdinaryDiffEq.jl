# Adding a new package to the common interface

DiffEq's distributed infrastructure allows anyone to add new packages to the
common interface. This set of the documentation explains how this is done.
An example package is [DASRK.jl](https://github.com/JuliaDiffEq/DASKR.jl/blob/master/src/common.jl)
whose full common interface bindings are contained in `common.jl`.

## Defining the types

You should start by defining a common supertype for all your algorithm types.
DASKR has DAE algorithms, so it defines

```julia
abstract type DASKRDAEAlgorithm{LinearSolver} <: DiffEqBase.AbstractDAEAlgorithm end
```

that its algorithms are all `AbstractDAEAlgorithm`s and gives them a possible
type parameter as well. Then the concrete algorithms are specified. Special
options (i.e. non-common interface options) for the algorithm go in this type.
Here there is a choice for a linear solver internally, so we allow the user
to set this:

```julia
struct daskr{LinearSolver} <: DASKRDAEAlgorithm{LinearSolver} end
daskr(; linear_solver = :Dense) = daskr{linear_solver}()
```

In many (most?) cases, no extra constructor is needed since there are no
extra options.

## Defining the overloads

Now you need to add DiffEqBase. The package should reexport `DiffEqBase.jl` to
make using it alone act naturally, and this is done by:

```julia
using Reexport
@reexport using DiffEqBase
```

now we overload `__solve` from `DiffEqBase.jl` to act on our algorithm. Here's a
possible signature:

```julia
function DiffEqBase.__solve{uType, duType, tType, isinplace, LinearSolver}(
        prob::AbstractDAEProblem{uType, duType, tType, isinplace},
        alg::DASKRDAEAlgorithm{LinearSolver},
        timeseries = [], ts = [], ks = []; verbose = true,
        callback = nothing, abstol = 1 / 10^6, reltol = 1 / 10^3,
        saveat = Float64[], adaptive = true, maxiters = Int(1e5),
        timeseries_errors = true, save_everystep = isempty(saveat), dense = save_everystep,
        save_start = true, save_timeseries = nothing,
        userdata = nothing,
        kwargs...)
    # do something
end
```

Basically you just dispatch on your algorithm supertype (and refine the `Problem`
choice as well so it errors on the wrong problem types), then list the common
interface options. `timeseries = [], ts = [], ks = [];` are for pre-allocation
of arrays, and this may change in the near future but you can use this to pre-allocate
the `t`, `u`, and `k` arrays (`k` being internal steps if saved).

## Defining the solution

In `solve` you do option handling and call your solver. At the end, you return
the solution via:

```julia
build_solution(prob, alg, ts, timeseries,
    du = dures,
    dense = dense,
    timeseries_errors = timeseries_errors,
    retcode = :Success)
```

Giving `du` is only currently allowed for DAEs and is optional. The errors
flags (`timeseries_errors` and `dense_errors`) are flags for allowing the solution
to calculate errors which should be passed through the `solve` function if
applicable. A proper `retcode` should be placed or it will default to `:Default`.

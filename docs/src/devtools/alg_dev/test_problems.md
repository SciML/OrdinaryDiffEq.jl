# Test Problems

For every problem, one can turn it into a test problem by adding `analytical`
to the `DEFunction`.

## No Analytical Solution

However, in many cases the analytical solution cannot be found, and therefore
one uses a low-tolerance calculation as a stand-in for a solution. The JuliaDiffEq
ecosystem supports this through the `TestSolution` type in DiffEqDevTools. There
are three constructors. The code is simple, so here it is:

```julia
mutable struct TestSolution <: DESolution
    t::Any
    u::Any
    interp::Any
    dense::Any
end
(T::TestSolution)(t) = T.interp(t)
TestSolution(t, u) = TestSolution(t, u, nothing, false)
TestSolution(t, u, interp) = TestSolution(t, u, interp, true)
TestSolution(interp::DESolution) = TestSolution(nothing, nothing, interp, true)
```

This acts like a solution. When used in conjunction with `apprxtrue`:

```julia
appxtrue(sol::AbstractODESolution, sol2::TestSolution)
```

you can use it to build a `TestSolution` from a problem (like `ODETestSolution`)
which holds the errors  If you only give it `t` and `u`, then it can only calculate
the final error. If the `TestSolution` has an interpolation, it will define
timeseries and dense errors.

(Note: I would like it so that way the timeseries
error will be calculated on the times of `sol.t in sol2.t` which would act nicely
with `tstops` and when interpolations don't exist, but haven't gotten to it!)

These can then be passed to other functionality. For example, the benchmarking
functions allow one to set `appxsol` which is a `TestSolution` for the benchmark
solution to calculate errors against, and `error_estimate` allows one to choose
which error estimate to use in the benchmarking (defaults to `:final`).

## Related Functions

```@docs
DiffEqDevTools.appxtrue
```

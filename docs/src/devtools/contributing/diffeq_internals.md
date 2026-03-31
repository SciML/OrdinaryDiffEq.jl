# The DiffEq Internals

The DiffEq solvers, OrdinaryDiffEq, StochasticDiffEq, FiniteElementDiffEq, etc.
all follow a similar scheme which leads to rapid development and high performance.
This portion of the documentation explains how the algorithms are written.

## Developing New Solver Algorithms

The easiest way to get started would be to add new solver algorithms. This is a
pretty simple task as there are tools which put you right into the "hot loop".
For example, take a look at the ODE solver code. The mode `solve(::ODEProblem,::OrdinaryDiffEqAlgorithm)`
is glue code to a bunch of solver algorithms. The algorithms which are coded
in DifferentialEquations.jl can be found in src/integrators.jl. The actual step
is denoted by the `perform_step!(integrator)` function. For example,
take a look at the Midpoint method's implementation:

```julia
@inline function perform_step!(integrator, cache::MidpointConstantCache, f = integrator.f)
    @unpack t, dt, uprev, u, k = integrator
    halfdt = dt / 2
    k = integrator.fsalfirst
    k = f(t + halfdt, uprev + halfdt * k)
    u = uprev + dt * k
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    @pack integrator = t, dt, u
end
```

The available items are all unloaded from the `integrator` in the first line.
`fsalfirst` inherits the value of `fsallast` on the last line. The algorithm is
written in this form so that way the derivative of both endpoints is defined, allowing
the vector `integrator.k` to determine a  Hermite interpolation polynomial (in general,
the `k` values for each algorithm form the interpolating polynomial). Other than that,
the algorithm is as basic as it gets for the Midpoint method, making sure to set
`fsallast` at the last line. The results are then packaged back into the integrator
to mutate the state.

Note the name `ConstantCache`. In OrdinaryDiffEq.jl, the `Algorithm` types are used
for holding type information about which solver to choose, and its the "cache" types
which then hold the internal caches and state. `ConstantCache` types are for
non in-place calculations `f(u,p,t)` and `Cache` types (like `MidpointCache`) include
all of the internal arrays. Their constructor is specified in the `cache.jl` file.
The cache (and the first `fsalfirst`) is initialized in the `initialize!` function
next to the cache's `perform_step!` function.

The main inner loop can be summarized by the `solve!` command:

```julia
@inbounds while !isempty(integrator.opts.tstops)
    while integrator.tdir * integrator.t < integrator.tdir * top(integrator.opts.tstops)
        loopheader!(integrator)
        @ode_exit_conditions
        perform_step!(integrator, integrator.cache)
        loopfooter!(integrator)
        if isempty(integrator.opts.tstops)
            break
        end
    end
    handle_tstop!(integrator)
end
postamble!(integrator)
```

The algorithm runs until `tstop` is empty. It hits the `loopheader!` in order to
accept/reject the previous step and choose a new `dt`. This is done at the top
so that way the iterator interface happens "mid-step" instead of "post-step".
In here the algorithm is also chosen for the switching algorithms.

Then the `perform_step!` is called. (The exit conditions throw an error if necessary.)
Afterwards, the `loopfooter!` is used to calculate new timesteps, save, and apply the
callbacks. If a value of `tstops` is hit, the algorithm breaks out of the inner-most
loop to save the value.

Adding algorithms to the other problems is very similar, just in a different package.

## Extras

If the method is a FSAL method then it needs to be set via `isfsal` and `fsalfirst`
should be defined before the loop, with `fsallast` what's pushed up to `fsalfirst`
upon a successful step. See `:DP5` for an example.

If tests fail due to units (i.e. Unitful), don't worry. I would be willing to fix
that up. To do so, you have to make sure you keep separate your `rateType`s and
your `uType`s since the rates from `f` will have units of `u` but divided by
a unit of time. If you simply try to write these into `u`, the units part will
fail (normally you have to multiply by a ``dt``).

If you want to access the value of u at the second-last time point, you can use uprev2 value. But in order to copy uprev value to uprev2 after each timestep, you need to make `alg_extrapolates(alg::Your_Alg) = true` in the alg_utils.jl file.

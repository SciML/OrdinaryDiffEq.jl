# Benchmark Suite

DifferentialEquations.jl provides a benchmarking suite to be able to test the
difference in error, speed, and efficiency between algorithms. DifferentialEquations.jl includes current benchmarking notebooks to help users
understand the performance of the methods. These benchmarking notebooks use the included benchmarking suite. There are two parts to the benchmarking suite: shootouts and work-precision. The `Shootout` tests methods head-to-head for timing and error on the same problem. A `WorkPrecision` draws a work-precision diagram
for the algorithms in question on the chosen problem.

### Rendered Benchmarks

The rendered SciML Benchmarks can be found at [benchmarks.sciml.ai](https://benchmarks.sciml.ai/stable/). The source code
for the benchmarks can be found at [https://github.com/SciML/SciMLBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl).

### Shootout

A
shootout is where you compare between algorithms. For example, to see how
different Runge-Kutta algorithms fair against each other, one can define a setup
which is a dictionary of Symbols to Any, where the symbol is the keyword argument.
Then you call `Shootout` on that setup. The code is as follows:

```julia
using OrdinaryDiffEq, DiffEqProblemLibrary.ODEProblemLibrary, DiffEqDevTools, ODE,
    ODEInterface, ODEInterfaceDiffEq

ODEProblemLibrary.importodeproblems()
prob = ODEProblemLibrary.prob_ode_2Dlinear
setups = [
    Dict(:alg => DP5())
    Dict(:abstol => 1.0e-3, :reltol => 1.0e-6, :alg => ode45()) # Fix ODE to be normal
    Dict(:alg => dopri5())
]
names = ["DifferentialEquations"; "ODE"; "ODEInterface"]
shoot = Shootout(prob, setups; dt = 1 / 2^(10), names = names)
```

Note that keyword arguments applied to `Shootout` are applied to every run, so
in this example every run has the same starting timestep.  Here we explicitly chose names.
If you don't, then the algorithm name is the default.
This returns a Shootout type which holds the times it took for each algorithm
and the errors. Using these, it calculates the efficiency defined as
1/(error*time), i.e. if the error is low or the run was quick then
it's efficient. `print(shoot)` will show all of this information,
and `plot(shoot)` will show the efficiencies of the algorithms
in comparison to each other.

For every benchmark function there is a special keyword `numruns` which controls
the number of runs used in the time estimate. To be more precise, these functions
by default run the algorithm 20 times on the problem and take the average time.
This amount can be increased and decreased as needed.

The keyword `appxsol` allows for specifying a reference against which the error is computed.
The method of error computation can be specified by the keyword `error_estimate` with values `:L2` for the L2 error over the solution time interval, `:l2` calculates the l2 error at the actual steps and the default `:final` only compares the endpoints.

A ShootoutSet is a where you define a vector of probs and tspans and run a shootout
on each of these values.

### WorkPrecision

A WorkPrecision calculates the necessary components of a work-precision plot. This
shows how time scales with the user chosen tolerances on a given problem. To make
a WorkPrecision, you give it a vector of absolute and relative tolerances:

```julia
abstols = 1 ./ 10 .^ (3:10)
reltols = 1 ./ 10 .^ (3:10)
wp = WorkPrecision(prob, DP5(), abstols, reltols; name = "Dormand-Prince 4/5")
```

If we want to plot many WorkPrecisions together in order to compare between
algorithms, you can make a WorkPrecisionSet. To do so, you pass the setups
into the function as well:

```julia
wp_set = WorkPrecisionSet(prob, tspan, abstols, reltols, setups; numruns = 2)
setups = [
    Dict(:alg => RK4()); Dict(:alg => Euler()); Dict(:alg => BS3());
    Dict(:alg => Midpoint()); Dict(:alg => BS5()); Dict(:alg => DP5())
]
wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)
```

Both of these types have a plot recipe to produce a work-precision diagram,
and a print which will show some relevant information.

### Tags and comparison plots

A benchmark usually wants several views of the same data: each family of methods on its
own, then the best of each family against each other, with a couple of reference methods
in every plot. Preset tags derived from algorithm traits and supertypes produce all of
those from a single run. Add only benchmark-specific tags such as `:reference`, request
every error metric the plots need with `error_estimates`, and slice the result afterwards:

```julia
setups = [
    Dict(:alg => Rosenbrock23()),
    Dict(:alg => Rodas5P()),
    Dict(:alg => TRBDF2()),
    Dict(:alg => KenCarp4()),
    Dict(:alg => RadauIIA5(), :tags => [:reference]),
]
wp_set = WorkPrecisionSet(
    prob, abstols, reltols, setups;
    error_estimates = [:final, :l2], appxsol = test_sol
)

plot(wp_set, tags = [:rosenbrock])                     # one family
plot(best_of_families(wp_set, [:rosenbrock, :sdirk]))  # cross-family comparison
plot(wp_set, x = :l2)                                  # a second error metric, no re-solve

# a family against the baseline
plot(
    wp_set, tags = [:sdirk], include_tags = [:reference],
    reference_tags = [:reference]
)
```

For example, `auto_tags(KenCarp4())` includes `:order_4`, `:adaptive`, `:implicit`,
`:sdirk`, `:esdirk`, and `:split`. Explicit setup tags are appended without duplicates.
Use `:auto_tags => false` when a setup needs only its manually supplied tags.
[`tag_kind`](@ref) distinguishes algorithm families from traits, benchmark roles,
providers, variants, and problem domains. This lets `autoplot(wp_set)` discover family
views without treating tags such as `:order_4` or `:reference` as families. Pass
`families` explicitly for one-off custom family tags.

`plot(wp_set; tags)` keeps the entries carrying all of `tags`, `include_tags` adds
entries back regardless of that filter, and `exclude_tags` drops entries. Entries
matching `reference_tags` are drawn in a separate, de-emphasized style controlled by
`reference_style`. [`autoplot`](@ref) returns the whole standard collection of subsets
at once, keyed by name.

Slow configurations can be capped with `timeout` (seconds per tolerance): a solve is
never interrupted, but once one exceeds the budget its repeated timing runs are skipped
and the point is recorded as `NaN`, which the plot recipe drops.

## API

```@docs
DiffEqDevTools.Shootout
DiffEqDevTools.ShootoutSet
DiffEqDevTools.WorkPrecision
DiffEqDevTools.WorkPrecisionSet
DiffEqDevTools.get_sample_errors
```

### Tagging and comparison helpers

```@docs
DiffEqDevTools.auto_tags
DiffEqDevTools.tag_kind
DiffEqDevTools.get_tags
DiffEqDevTools.unique_tags
DiffEqDevTools.filter_by_tags
DiffEqDevTools.exclude_by_tags
DiffEqDevTools.merge_wp_sets
DiffEqDevTools.available_errors
DiffEqDevTools.wp_area
DiffEqDevTools.best_by_tag
DiffEqDevTools.best_of_families
DiffEqDevTools.with_autodiff_variants
DiffEqDevTools.autoplot
```

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
setups = [Dict(:alg => DP5())
          Dict(:abstol => 1e-3, :reltol => 1e-6, :alg => ode45()) # Fix ODE to be normal
          Dict(:alg => dopri5())]
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

A WorkPrecision calculates the necessary componnets of a work-precision plot. This
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
setups = [Dict(:alg => RK4()); Dict(:alg => Euler()); Dict(:alg => BS3());
          Dict(:alg => Midpoint()); Dict(:alg => BS5()); Dict(:alg => DP5())]
wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)
```

Both of these types have a plot recipe to produce a work-precision diagram,
and a print which will show some relevant information.

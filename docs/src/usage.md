# Usage

OrdinaryDiffEq.jl is part of the SciML common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes an OrdinaryDiffEq.jl algorithm to `solve`. For example, we can solve the [ODE tutorial from the docs](https://docs.sciml.ai/DiffEqDocs/dev/getting_started/) using the `Tsit5()` algorithm:

```julia
using OrdinaryDiffEq
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
using Plots
plot(
    sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!", # legend = false
)
plot!(sol.t, t -> 0.5 * exp(1.01 * t), lw = 3, ls = :dash, label = "True Solution!")
```
`Tsit5()` is a good default choice for many non-stiff ODEs. For stiff problems,
consider using an implicit method such as `Rodas5P()` or `TRBDF2()`.

That example uses the out-of-place syntax `f(u,p,t)`, while the in-place syntax (more efficient for systems of equations) is shown in the Lorenz example:

```julia
using OrdinaryDiffEq
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())
using Plots;
plot(sol, idxs = (1, 2, 3));
```

Very fast static array versions can be specifically compiled to the size of your model. For example:

```julia
using OrdinaryDiffEq, StaticArrays
function lorenz(u, p, t)
    return SA[10.0 * (u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end
u0 = SA[1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, Tsit5())
```

For “refined ODEs”, like dynamical equations and `SecondOrderODEProblem`s, refer to the [DiffEqDocs](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/). For example, in [DiffEqTutorials.jl](https://github.com/SciML/SciMLTutorials.jl) we show how to solve equations of motion using symplectic methods:

```julia
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dv[1] = dx = -x - 2 * x * y
    dv[2] = dy = y^2 - y - x^2
    return
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol2 = solve(prob, KahanLi8(), dt = 1 / 10);
```

Other refined forms are IMEX and semi-linear ODEs (for exponential integrators).

## Reactant compilation

An explicit ODE solve can be part of a [`Reactant.@jit`](https://enzymead.github.io/Reactant.jl/stable/api/#Reactant.@jit) compiled function. Both adaptive and fixed-step solver loops are staged as device-side loops, so the compiled executable can be reused with new state and parameter values:

```julia
using OrdinaryDiffEq, Reactant

f(u, p, t) = p .* u

function compiled_solve(u, p)
    prob = ODEProblem(f, u, (0.0f0, 1.0f0), p)
    return solve(prob, Tsit5())
end

u0 = Reactant.to_rarray(Float32[1, 2])
p = Reactant.to_rarray(Float32[-1])
sol = Reactant.@jit compiled_solve(u0, p)
Array(sol.u[end])
```

Reactant requires statically shaped outputs. A compiled solve therefore returns an endpoint-only `ODESolution`: `sol.u` and `sol.t` contain the final state and time, while `sol.prob`, `sol.stats`, and `sol.interp` are `nothing`. Adaptive solves currently require a `PIController`; the tested algorithms are the `Tsit5` and Verner explicit Runge–Kutta families. Implicit algorithms, saving intermediate or partial states (`saveat` or `save_idxs`), callbacks, user `tstops`, discontinuity handling, `force_dtmin`, progress reporting, custom domain or instability checks, and step limiters are not currently supported inside Reactant compilation and produce an `ArgumentError` instead of silently changing the solve.

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl ODE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/), [Dynamical ODE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dynamical_solve/), and the [Split ODE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/split_ode_solve/) pages.

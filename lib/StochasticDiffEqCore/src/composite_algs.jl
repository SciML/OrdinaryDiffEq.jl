mutable struct AutoSwitchCache{nAlg, sAlg, tolType, T}
    count::Int
    successive_switches::Int
    nonstiffalg::nAlg
    stiffalg::sAlg
    is_stiffalg::Bool
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
end

"""
    AutoSwitch(nonstiffalg, stiffalg; maxstiffstep = 10, maxnonstiffstep = 3,
               nonstifftol = 9//10, stifftol = 9//10, dtfac = 2,
               stiffalgfirst = false, switch_max = 5)

Choice function that switches a [`StochasticCompositeAlgorithm`](@ref) between a
nonstiff and a stiff method based on an online stiffness estimate.

Each step the estimate `abs(eigen_est * dt / alg_stability_size(alg))` is compared
against a tolerance. Consecutive verdicts are counted, and the algorithm is only
switched once the count exceeds the corresponding threshold, which prevents
thrashing on a borderline problem.

## Keyword Arguments

  - `maxstiffstep`: Number of consecutive stiff verdicts required before switching to
    `stiffalg` (default: `10`).
  - `maxnonstiffstep`: Number of consecutive nonstiff verdicts required before
    switching back to `nonstiffalg` (default: `3`).
  - `nonstifftol`, `stifftol`: Stiffness-ratio tolerances used while the nonstiff and
    the stiff algorithm is active, respectively (default: `9//10` for both).
  - `dtfac`: Factor applied to `dt` at a switch — `dt` is multiplied by it when moving
    to the stiff method and divided by it when moving back (default: `2`).
  - `stiffalgfirst`: Start with `stiffalg` instead of `nonstiffalg` (default: `false`).
  - `switch_max`: Number of successive nonstiff verdicts after which the error check
    is re-enabled (default: `5`).

Use [`AutoAlgSwitch`](@ref) to build the composite algorithm directly.
"""
struct AutoSwitch{nAlg, sAlg, tolType, T}
    nonstiffalg::nAlg
    stiffalg::sAlg
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
end
function AutoSwitch(
        nonstiffalg, stiffalg; maxstiffstep = 10, maxnonstiffstep = 3,
        nonstifftol = 9 // 10, stifftol = 9 // 10, dtfac = 2, stiffalgfirst = false,
        switch_max = 5
    )
    return AutoSwitch(
        nonstiffalg, stiffalg, maxstiffstep, maxnonstiffstep,
        promote(nonstifftol, stifftol)..., dtfac, stiffalgfirst, switch_max
    )
end

function is_stiff(integrator, alg, ntol, stol, is_stiffalg)
    eigen_est, dt = integrator.eigen_est, integrator.dt
    stiffness = abs(eigen_est * dt / alg_stability_size(alg)) # `abs` here is just for safety
    tol = is_stiffalg ? stol : ntol
    os = oneunit(stiffness)
    bool = stiffness > os * tol

    if !bool
        integrator.alg.choice_function.successive_switches += 1
    else
        integrator.alg.choice_function.successive_switches = 0
    end

    integrator.do_error_check = integrator.alg.choice_function.successive_switches >
        integrator.alg.choice_function.switch_max || !bool
    return bool
end

function (AS::AutoSwitchCache)(integrator)
    integrator.iter == 0 && return Int(AS.stiffalgfirst) + 1
    dt = integrator.dt
    # Successive stiffness test positives are counted by a positive integer,
    # and successive stiffness test negatives are counted by a negative integer
    AS.count = is_stiff(
            integrator, AS.nonstiffalg, AS.nonstifftol, AS.stifftol, AS.is_stiffalg
        ) ?
        AS.count < 0 ? 1 : AS.count + 1 :
        AS.count > 0 ? -1 : AS.count - 1
    if (!AS.is_stiffalg && AS.count > AS.maxstiffstep)
        integrator.dt = dt * AS.dtfac
        AS.is_stiffalg = true
    elseif (AS.is_stiffalg && AS.count < -AS.maxnonstiffstep)
        integrator.dt = dt / AS.dtfac
        AS.is_stiffalg = false
    end
    return Int(AS.is_stiffalg) + 1
end

"""
    AutoAlgSwitch(nonstiffalg, stiffalg; kwargs...) -> StochasticCompositeAlgorithm

Build a stiffness-switching composite of `nonstiffalg` and `stiffalg`.

This is the constructor behind the `Auto*` SDE solvers (for example `AutoSOSRI2`):
it wraps the two algorithms in a [`StochasticCompositeAlgorithm`](@ref) whose choice
function is an [`AutoSwitch`](@ref). All keyword arguments are forwarded to
`AutoSwitch`.

```julia
alg = AutoAlgSwitch(SOSRI(), ImplicitEM(); maxstiffstep = 5)
```
"""
function AutoAlgSwitch(nonstiffalg, stiffalg; kwargs...)
    AS = AutoSwitch(nonstiffalg, stiffalg; kwargs...)
    return StochasticCompositeAlgorithm((nonstiffalg, stiffalg), AS)
end

# AutoSOSRI2/AutoSOSRA2 are defined in the StochasticDiffEq umbrella package
# since they reference concrete types from multiple solver subpackages.

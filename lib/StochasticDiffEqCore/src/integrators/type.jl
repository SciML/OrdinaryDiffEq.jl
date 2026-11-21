"""
    SDEAlgTypes

Union of the algorithm supertypes that the SDE integrator loop handles, i.e.
`StochasticDiffEqAlgorithm` and `StochasticDiffEqRODEAlgorithm`.

Used to constrain [`SDEIntegrator`](@ref) and to write methods that apply to both SDE
and RODE algorithms at once.
"""
const SDEAlgTypes = Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}

"""
    SDEIntegrator

The integrator type for SDE and RODE problems: an `ODEIntegrator` specialized to an
algorithm in [`SDEAlgTypes`](@ref).

SDE integration reuses OrdinaryDiffEqCore's integrator machinery and adds the noise
process in the `W` field (and `P` for jump problems), so this alias is what
SDE-specific methods — resizing, noise cache handling, the interpolation and step
interface — dispatch on.
"""
const SDEIntegrator = ODEIntegrator{<:SDEAlgTypes}

"""
    SDEOptions

Alias for `OrdinaryDiffEqCore.DEOptions`, the container holding the solver options
(tolerances, callbacks, save settings) of a running integrator.

Kept as its own name for downstream packages such as StochasticDelayDiffEq that were
written against the pre-migration `StochasticDiffEq.SDEOptions`.
"""
const SDEOptions = OrdinaryDiffEqCore.DEOptions

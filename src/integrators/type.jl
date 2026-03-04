const SDEAlgTypes = Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
const SDEIntegrator = ODEIntegrator{<:SDEAlgTypes}
# Backwards-compatibility alias for downstream packages (e.g. StochasticDelayDiffEq)
const SDEOptions = OrdinaryDiffEqCore.DEOptions

const SDEAlgTypes = Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
const SDEIntegrator = ODEIntegrator{<:SDEAlgTypes}

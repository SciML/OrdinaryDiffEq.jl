# the following are setup per how integrators are implemented in OrdinaryDiffEq and
# StochasticDiffEq and provide dispatch points that JumpProcesses and others can use.

"""
    get_tstops(integrator) -> Any

Return the timestep-stop data structure owned by `integrator`.

Integrator implementations specialize this accessor so callback and jump-process
machinery can inspect pending `tstops` without depending on integrator fields.
"""
function get_tstops(integ::DEIntegrator)
    error("get_tstops not implemented for integrators of type $(nameof(typeof(integ)))")
end

"""
    get_tstops_array(integrator) -> AbstractVector

Return the array-like storage containing pending timestep stops for `integrator`.
"""
function get_tstops_array(integ::DEIntegrator)
    error("get_tstops_array not implemented for integrators of type $(nameof(typeof(integ)))")
end

"""
    get_tstops_max(integrator)

Return the largest pending timestep stop for `integrator`.
"""
function get_tstops_max(integ::DEIntegrator)
    error("get_tstops_max not implemented for integrators of type $(nameof(typeof(integ)))")
end

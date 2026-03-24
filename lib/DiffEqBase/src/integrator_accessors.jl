# the following are setup per how integrators are implemented in OrdinaryDiffEq and
# StochasticDiffEq and provide dispatch points that JumpProcesses and others can use.

function get_tstops(integ::DEIntegrator)
    error("get_tstops not implemented for integrators of type $(nameof(typeof(integ)))")
end
function get_tstops_array(integ::DEIntegrator)
    error("get_tstops_array not implemented for integrators of type $(nameof(typeof(integ)))")
end
function get_tstops_max(integ::DEIntegrator)
    error("get_tstops_max not implemented for integrators of type $(nameof(typeof(integ)))")
end

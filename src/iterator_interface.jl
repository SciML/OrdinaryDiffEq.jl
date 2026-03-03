# Override step! for SDE integrators to return nothing (matching old behavior).
# ODE's _step! returns integrator.sol.retcode, but on Julia 1.10 the ReturnCode.T
# return value gets boxed (16 bytes), causing allocation test failures.
# SDE's old step! always returned nothing, so we preserve that behavior.
@inline function DiffEqBase.step!(integrator::SDEIntegrator)
    _step!(integrator)
    return nothing
end

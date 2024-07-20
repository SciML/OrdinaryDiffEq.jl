# checks whether the controller should accept a step based on the error estimate
@inline function accept_step_controller(integrator, controller::AbstractController)
    return integrator.EEst <= 1
end

@inline function accept_step_controller(integrator, controller::PIDController)
    return integrator.qold >= controller.accept_safety
end

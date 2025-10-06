function resize_non_user_cache!(integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i)
    resize_J_W!(cache, integrator, i)
    resize_jac_config!(cache, integrator)
    resize_grad_config!(cache, integrator)

    nothing
end

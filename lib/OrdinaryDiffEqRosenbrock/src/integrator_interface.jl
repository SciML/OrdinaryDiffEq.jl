function resize_non_user_cache!(
        integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i
    )
    resize_J_W!(cache, integrator, i, cache.alg)
    resize_jac_config!(cache, integrator, cache.alg)
    resize_grad_config!(cache, integrator, cache.alg)

    return nothing
end

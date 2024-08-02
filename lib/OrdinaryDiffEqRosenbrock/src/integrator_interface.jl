function resize_non_user_cache!(integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)
    resize_jac_config!(cache.jac_config, i)
    resize_grad_config!(cache.grad_config, i)
    nothing
end
function resize_non_user_cache!(integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)

    uf = cache.uf
    uf = SciMLBase.@set uf.f = SciMLBase.unwrapped_f(uf.f)

    cache.jac_config = resize_jac_config!(uf, cache.du1, cache.jac_config, alg_autodiff(integrator.alg), integrator.u)

    if alg_autodiff(integrator.alg) isa AutoSparse
        ad = ADTypes.dense_ad(alg_autodiff(integrator.alg))
    else
        ad = alg_autodiff(integrator.alg)
    end

    cache.grad_config = resize_grad_config!(cache.tf, cache.du1, cache.grad_config, ad, integrator.t)

    nothing
end

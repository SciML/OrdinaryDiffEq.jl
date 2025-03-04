function resize_non_user_cache!(integrator::ODEIntegrator,
        cache::RosenbrockMutableCache, i)
    cache.J = similar(cache.J, i, i)
    cache.W = similar(cache.W, i, i)

    uf = cache.uf
    uf = SciMLBase.@set uf.f = SciMLBase.unwrapped_f(uf.f)

    if !isnothing(cache.grad_config) && !isnothing(cache.grad_config[1])

        if alg_autodiff(integrator.alg) isa AutoSparse
            ad = ADTypes.dense_ad(alg_autodiff(integrator.alg))
        else
            ad = alg_autodiff(integrator.alg)
        end

        if ad isa AutoFiniteDiff
            ad_right = SciMLBase.@set ad.dir = 1
            ad_left = SciMLBase.@set ad.dir = -1
        else
            ad_right = ad
            ad_left = ad
        end

        cache.grad_config = ([resize_grad_config!(cache.tf, cache.du1, config, ad, integrator.t) for (ad, config) in zip((ad_right, ad_left), cache.grad_config)]...,)
    end

    if !isnothing(cache.jac_config) && !isnothing(cache.jac_config[1])

        # for correct FiniteDiff dirs
        autodiff_alg = alg_autodiff(integrator.alg)
        if autodiff_alg isa AutoFiniteDiff
            ad_right = SciMLBase.@set autodiff_alg.dir = 1
            ad_left = SciMLBase.@set autodiff_alg.dir = -1
        else
            ad_right = autodiff_alg
            ad_left = autodiff_alg
        end

        cache.jac_config = ([resize_jac_config!(
                                   uf, cache.du1, config, ad, integrator.u)
                               for (ad, config) in zip(
            (ad_right, ad_left), cache.jac_config)]...,)
    end

    nothing
end

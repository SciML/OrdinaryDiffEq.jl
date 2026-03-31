#methods from https://doi.org/10.1016/j.cam.2009.11.010
@muladd function perform_step!(integrator, cache::WangLi3SMil_AConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f(uprev, p, t)
    u = uprev + dt * k

    #stage 2
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / (integrator.sqdt)
    u = u - (dt / 2) * k1

    #stage 3
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_ACache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f(k, uprev, p, t)
    @.. u = uprev + dt * k

    #stage 2
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / (integrator.sqdt)
    @.. u = u - (dt / 2) * k1

    #stage 3
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u + k * W.dW + (W.dW^2 / 2) * k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_BConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f.g(uprev, p, t)
    tmp = uprev + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / integrator.sqdt
    u = uprev .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1

    #stage 2
    k = integrator.f(u, p, t)
    u = u + dt * k

    #stage 3
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u - (dt / 2) * k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_BCache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f.g(k, uprev, p, t)
    @.. tmp = uprev + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / integrator.sqdt
    @.. u = uprev + W.dW * k + (W.dW^2 / 2) * k1

    #stage 2
    integrator.f(k, u, p, t)
    @.. u = u + dt * k

    #stage 3
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u - (dt / 2) * k1

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_CConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f(uprev, p, t)
    u = uprev + dt * k

    #stage 2
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / (integrator.sqdt)
    u = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1

    #stage 3
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u - (dt / 2) * k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_CCache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f(k, uprev, p, t)
    @.. u = uprev + dt * k

    #stage 2
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / (integrator.sqdt)
    @.. u = u + k * W.dW + (W.dW^2 / 2) * k1

    #stage 3
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u - (dt / 2) * k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_DConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f.g(uprev, p, t)
    tmp = uprev + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / (integrator.sqdt)
    u = uprev - (dt / 2) * k1

    #stage 2
    k = integrator.f(u, p, t)
    u = u + dt * k

    #stage 3
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_DCache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f.g(k, uprev, p, t)
    @.. tmp = uprev + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / (integrator.sqdt)
    @.. u = uprev - (dt / 2) * k1

    #stage 2
    integrator.f(k, u, p, t)
    @.. u = u + dt * k

    #stage 3
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u + k * W.dW + (W.dW^2 / 2) * k1
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_EConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f.g(uprev, p, t)
    tmp = uprev + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / (integrator.sqdt)
    u = uprev - (dt / 2) * k1

    #stage 2
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1

    #stage 3
    k = integrator.f(u, p, t)
    u = u + dt * k
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_ECache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f.g(k, uprev, p, t)
    @.. tmp = uprev + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / (integrator.sqdt)
    @.. u = uprev - (dt / 2) * k1

    #stage 2
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u + k * W.dW + (W.dW^2 / 2) * k1

    #stage 3
    integrator.f(k, u, p, t)
    @.. u = u + dt * k
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_FConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    #stage 1
    k = integrator.f.g(uprev, p, t)
    tmp = uprev + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / integrator.sqdt
    u = uprev .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k1

    #stage 2
    k = integrator.f.g(u, p, t)
    tmp = u + k * integrator.sqdt
    k1 = integrator.f.g(tmp, p, t)
    k1 = (k1 - k) / abs(integrator.sqdt)
    u = u - (dt / 2) * k1

    #stage 3
    k = integrator.f(u, p, t)
    u = u + dt * k
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::WangLi3SMil_FCache)
    (; tmp, k, k1) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    #stage 1
    integrator.f.g(k, uprev, p, t)
    @.. tmp = uprev + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / integrator.sqdt
    @.. u = uprev + W.dW * k + (W.dW^2 / 2) * k1

    #stage 2
    integrator.f.g(k, u, p, t)
    @.. tmp = u + k * integrator.sqdt
    integrator.f.g(k1, tmp, p, t)
    @.. k1 = (k1 - k) / abs(integrator.sqdt)
    @.. u = u - (dt / 2) * k1

    #stage 3
    integrator.f(k, u, p, t)
    @.. u = u + dt * k
    integrator.u = u
end

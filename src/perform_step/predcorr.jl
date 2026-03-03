@muladd function perform_step!(integrator, cache::PCEulerConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    g = integrator.f.g
    (; theta, eta, ggprime) = integrator.alg
    dW = W.dW
    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        gdW_n = g(uprev, p, t) * dW
    else
        gdW_n = g(uprev, p, t) .* dW
    end

    f_n = f(uprev, p, t)
    ubar = uprev + dt * f_n + gdW_n

    fbar_n = f_n - eta * ggprime(uprev, p, t)

    tnp1 = t + dt
    fbar_np1 = f(ubar, p, tnp1) - eta * ggprime(ubar, p, tnp1)

    if !is_diagonal_noise(integrator.sol.prob) || dW isa Number
        gdW_np1 = g(ubar, p, tnp1) * dW
    else
        gdW_np1 = g(ubar, p, tnp1) .* dW
    end

    u = uprev + (theta * dt) * fbar_np1 + ((1 - theta) * dt) * fbar_n + eta * gdW_np1 + (1 - eta) * gdW_n
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::PCEulerCache)
    (; utmp, ftmp, gtmp, gdWtmp, bbprimetmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    g = integrator.f.g
    (; theta, eta, ggprime) = integrator.alg
    dW = W.dW

    f(ftmp, uprev, p, t)
    g(gtmp, uprev, p, t)
    ggprime(bbprimetmp, uprev, p, t)

    if is_diagonal_noise(integrator.sol.prob)
        @.. gtmp *= dW
    else
        mul!(gdWtmp, gtmp, dW)
    end

    @.. utmp = uprev + ftmp * dt + gdWtmp
    @.. u = uprev + (1 - theta) * (ftmp - eta * bbprimetmp) * dt + (1 - eta) * gdWtmp

    tnp1 = t + dt
    f(ftmp, utmp, p, tnp1)
    g(gtmp, utmp, p, tnp1)
    ggprime(bbprimetmp, utmp, p, tnp1)

    if is_diagonal_noise(integrator.sol.prob)
        @.. gtmp *= dW
    else
        mul!(gdWtmp, gtmp, dW)
    end

    @.. u += theta * (ftmp - eta * bbprimetmp) * dt + eta * gdWtmp
end

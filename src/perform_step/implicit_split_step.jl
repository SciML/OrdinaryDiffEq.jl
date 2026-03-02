@muladd function perform_step!(
        integrator, cache::Union{
            ISSEMConstantCache, ISSEulerHeunConstantCache,
        }
    )
    (; t, dt, uprev, u, p, f) = integrator
    (; nlsolver) = cache
    alg = unwrap_alg(integrator, true)
    theta = alg.theta
    alg.symplectic ? a = dt / 2 : a = theta * dt
    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    # TODO: Stochastic extrapolants?
    u = uprev

    repeat_step = false

    L = integrator.f.g(uprev, p, t)
    ftmp = integrator.f(uprev, p, t)

    if alg.symplectic
        z = zero(u) # constant extrapolation, justified by ODE IM
    else
        z = dt * ftmp # linear extrapolation
    end
    nlsolver.z = z

    # nlsolver.c should be the Butcher tableau coefficient (time fraction), not coefficient*dt
    # OrdinaryDiffEqNonlinearSolve computes tstep = t + c * dt, so c should be in [0,1]
    nlsolver.c = alg.symplectic ? one(t) / 2 : theta
    if alg.symplectic
        # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
        #u = uprev + z/2
        tmp = uprev
    else
        tmp = uprev + dt * (1 - theta) * ftmp
    end
    nlsolver.tmp = tmp

    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    if alg.symplectic
        u = tmp + z
    else
        u = tmp + theta * z
    end

    if !is_diagonal_noise(integrator.sol.prob)
        gtmp = L * integrator.W.dW
    else
        gtmp = L .* integrator.W.dW
    end

    if cache isa ISSEulerHeunConstantCache
        utilde = u + gtmp
        if !is_diagonal_noise(integrator.sol.prob)
            gtmp = ((integrator.f.g(utilde, p, t) + L) / 2) * integrator.W.dW
        else
            gtmp = ((integrator.f.g(utilde, p, t) + L) / 2) .* integrator.W.dW
        end
    end

    u += gtmp

    if integrator.opts.adaptive
        if has_Wfact(f)
            # This means the Jacobian was never computed!
            J = f.jac(uprev, p, t)
        else
            J = OrdinaryDiffEqDifferentiation.calc_J(integrator, nlsolver.cache)
        end
        Ed = dt * (dt * J * ftmp) / 2

        if cache isa ISSEMConstantCache
            K = @.. uprev + dt * ftmp
            if !is_diagonal_noise(integrator.sol.prob)
                g_sized = norm(L, 2)
                utilde = @.. K + integrator.sqdt * g_sized
                gtmp2 = integrator.f.g(utilde, p, t)
                g_sized2 = norm(gtmp2, 2)
                ggprime = (g_sized2 - g_sized) / (integrator.sqdt)
                dW_cache = integrator.W.dW .^ 2 .- dt
                diff_tmp = integrator.opts.internalnorm(dW_cache, t)
                En = ggprime * diff_tmp / 2
            else
                utilde = @.. K + integrator.sqdt * L
                ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
                En = ggprime .* (integrator.W.dW .^ 2 .- dt) ./ 2
            end
        elseif cache isa ISSEulerHeunConstantCache
            if !is_diagonal_noise(integrator.sol.prob)
                g_sized = norm(L, 2)
                utilde = @.. uprev + g_sized * integrator.sqdt
                gtmp2 = integrator.f.g(utilde, p, t)
                g_sized2 = norm(gtmp2, 2)
                ggprime = (g_sized2 - g_sized) / (integrator.sqdt)
                dW_cache = integrator.W.dW .^ 2
                diff_tmp = integrator.opts.internalnorm(dW_cache, t)
                En = ggprime * diff_tmp / 2
            else
                utilde = @.. uprev + L * integrator.sqdt
                ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
                En = ggprime .* (integrator.W.dW .^ 2) ./ 2
            end
        end

        resids = calculate_residuals(
            Ed, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(resids, t)
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Union{ISSEMCache, ISSEulerHeunCache})
    (; t, dt, uprev, u, p, f) = integrator
    (; gtmp, gtmp2, dW_cache, nlsolver, k, dz) = cache
    (; z, tmp) = nlsolver

    J = (OrdinaryDiffEqCore.isnewton(nlsolver) ? nlsolver.cache.J : nothing)
    alg = unwrap_alg(integrator, true)
    alg.symplectic ? a = dt / 2 : a = alg.theta * dt
    dW = integrator.W.dW
    mass_matrix = integrator.f.mass_matrix
    theta = alg.theta
    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)

    repeat_step = false

    if integrator.success_iter > 0 && !integrator.u_modified &&
            alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
    elseif alg.extrapolant == :linear
        @.. u = uprev + integrator.fsalfirst * dt
    else # :constant
        copyto!(u, uprev)
    end

    integrator.f(tmp, uprev, p, t)
    integrator.f.g(gtmp, uprev, p, t)

    if alg.symplectic
        @.. z = zero(eltype(u)) # Justified by ODE solvers, constraint extrapolation when IM
    else
        @.. z = dt * tmp # linear extrapolation
    end

    # Handle noise computations

    if is_diagonal_noise(integrator.sol.prob)
        @.. gtmp2 = gtmp * dW
    else
        mul!(gtmp2, gtmp, dW)
    end

    ###
    # adaptivity part
    if integrator.opts.adaptive
        if has_Wfact(f)
            # This means the Jacobian was never computed!
            f.jac(J, uprev, p, t)
        else
            OrdinaryDiffEqDifferentiation.calc_J!(J, integrator, nlsolver.cache)
        end

        mul!(vec(z), J, vec(tmp))
        @.. k = dt * dt * z / 2
        # k is Ed
        # dz is En

        if !is_diagonal_noise(integrator.sol.prob)
            g_sized = norm(gtmp, 2)
        else
            g_sized = gtmp
        end
        # z is utilde above
        if cache isa ISSEMCache
            @.. z = uprev + dt * tmp + integrator.sqdt * g_sized
        elseif cache isa ISSEulerHeunCache
            @.. z = uprev + integrator.sqdt * g_sized
        end

        if cache isa ISSEMCache
            if !is_diagonal_noise(integrator.sol.prob)
                integrator.f.g(gtmp, z, p, t)
                g_sized2 = norm(gtmp, 2)
                @.. dW_cache = dW .^ 2 - dt
                diff_tmp = integrator.opts.internalnorm(dW_cache, t)
                En = (g_sized2 - g_sized) / (2integrator.sqdt) * diff_tmp
                @.. dz = En
            else
                integrator.f.g(gtmp2, z, p, t)
                g_sized2 = gtmp2
                @.. dz = (g_sized2 - g_sized) / (2integrator.sqdt) * (dW .^ 2 - dt)
            end
        elseif cache isa ISSEulerHeunCache
            if !is_diagonal_noise(integrator.sol.prob)
                integrator.f.g(gtmp, z, p, t)
                g_sized2 = norm(gtmp, 2)
                @.. dW_cache = dW .^ 2
                diff_tmp = integrator.opts.internalnorm(dW_cache, t)
                En = (g_sized2 - g_sized) / (2integrator.sqdt) * diff_tmp
                @.. dz = En
            else
                integrator.f.g(gtmp2, z, p, t)
                g_sized2 = gtmp2
                @.. dz = (g_sized2 - g_sized) / (2integrator.sqdt) * (dW .^ 2)
            end
        end
    end
    ###

    if alg.symplectic
        #@.. u = uprev + z/2
        @.. tmp = uprev
    else
        #@.. u = uprev + dt*(1-theta)*tmp + theta*z
        @.. tmp = uprev + dt * (1 - theta) * tmp
    end
    # nlsolver.c should be the Butcher tableau coefficient (time fraction), not coefficient*dt
    # OrdinaryDiffEqNonlinearSolve computes tstep = t + c * dt, so c should be in [0,1]
    nlsolver.c = alg.symplectic ? one(t) / 2 : theta
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return

    if alg.symplectic
        @.. u = uprev + z
    else
        #@.. u = uprev + dt*(1-theta)*tmp + theta*z
        @.. u = tmp + theta * z
    end

    ##############################################################################
    if cache isa ISSEulerHeunCache
        gtmp3 = cache.gtmp3
        @.. z = u + gtmp2
        integrator.f.g(gtmp3, z, p, t)
        @.. gtmp = (gtmp3 + gtmp) / 2
        if is_diagonal_noise(integrator.sol.prob)
            @.. gtmp2 = gtmp * dW
        else
            mul!(gtmp2, gtmp, dW)
        end
    end
    @.. u += gtmp2

    ##############################################################################
    if integrator.opts.adaptive
        calculate_residuals!(
            tmp, k, dz, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
end

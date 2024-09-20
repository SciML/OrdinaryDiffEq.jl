################################################################################

#### ROS2 type method

@ROS2(:init)
@ROS2(:performstep)

################################################################################

#### ROS23 type method

@ROS23(:init)
@ROS23(:performstep)

################################################################################

#### ROS34PW type method

@ROS34PW(:init)
@ROS34PW(:performstep)

################################################################################

#### ROS4 type method

@Rosenbrock4(:performstep)

#### Arbirary order method
function initialize!(integrator, cache::RosenbrockCombinedConstantCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    for i in 1:integrator.kshortsize
        integrator.k[i] = zero(integrator.u)
    end
end

@muladd function perform_step!(integrator, cache::RosenbrockCombinedConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;tf, uf) = cache
    (;A, C, gamma, c, d, H) = cache.tab

    # Precalculations
    dtC = C ./ dt
    dtd = dt .* d
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    # Initialize ks
    num_stages = size(A, 1)
    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    linsolve_tmp = @.. du + dtd[1] * dT
    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    # constant number for type stability make sure this is greater than num_stages
    ks = ntuple(Returns(k1), 10)
    # Loop for stages
    for stage in 2:num_stages
        u = uprev
        for i in 1:(stage - 1)
            u = @.. u + A[stage, i] * ks[i]
        end

        du = f(u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        # Compute linsolve_tmp for current stage
        linsolve_tmp = zero(du)
        if mass_matrix === I
            for i in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, i] * ks[i]
            end
        else
            for i in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, i] * ks[i]
            end
            linsolve_tmp = mass_matrix * linsolve_tmp
        end
        linsolve_tmp = @.. du + dtd[stage] * dT + linsolve_tmp

        ks = Base.setindex(ks, _reshape(W \ -_vec(linsolve_tmp), axes(uprev)), stage)
        integrator.stats.nsolve += 1
    end
    #@show ks
    u = u .+ ks[num_stages]

    if integrator.opts.adaptive
        atmp = calculate_residuals(ks[num_stages], uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] = zero(integrator.k[1])
        end
        for i in 1:num_stages
            for j in eachindex(integrator.k)
                integrator.k[j] = @.. integrator.k[j] + H[j, i] * ks[i]
            end
        end
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
            (integrator.EEst < 1.0)
             k2 = 0.5 * (uprev + u +
                   0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
             du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
             du = f(k2, p, t + dt / 2)
             OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
             if mass_matrix === I
                 du2 = du1 - du
             else
                 du2 = mass_matrix * du1 - du
             end
             EEst = norm(du2) / norm(integrator.opts.abstol .+ integrator.opts.reltol .* k2)
             integrator.EEst = max(EEst, integrator.EEst)
         end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::RosenbrockCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    resize!(integrator.k, integrator.kshortsize)
    for i in 1:integrator.kshortsize
        integrator.k[i] = cache.dense[i]
    end
end

@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, dT, J, W, uf, tf, ks, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; A, C, gamma, c, d, H) = cache.tab

    # Assignments
    sizeu = size(u)
    uidx = eachindex(integrator.uprev)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtC = C .* inv(dt)
    dtd = dt .* d
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    for stage in 2:length(ks)
        u .= uprev
        for i in 1:(stage - 1)
            @.. u += A[stage, i] * ks[i]
        end

        stage_limiter!(u, integrator, p, t + c[stage] * dt)
        f(du, u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        du1 .= 0
        if mass_matrix === I
            for i in 1:(stage - 1)
                @.. du1 += dtC[stage, i] * ks[i]
            end
        else
            for i in 1:(stage - 1)
                @.. du1 += dtC[stage, i] * ks[i]
            end
            mul!(_vec(du2), mass_matrix, _vec(du1))
            du1 .= du2
        end
        @.. linsolve_tmp = du + dtd[stage] * dT + du1

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        @.. $(_vec(ks[stage])) = -linres.u
        integrator.stats.nsolve += 1
    end
    du .= ks[end]
    u .+= ks[end]

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            @.. du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                    1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                    -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                    0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        calculate_residuals!(atmp, ks[end], uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] .= 0
        end
        for i in eachindex(ks)
            for j in eachindex(integrator.k)
                @.. integrator.k[j] += H[j, i] * ks[i]
            end
        end
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
            (integrator.EEst < 1.0)
             ks[2] = 0.5 * (uprev + u +
                   0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
             du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
             f(du, ks[2], p, t + dt / 2)
             OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
             if mass_matrix === I
                 @.. du2 = du1 - du
             else
                 mul!(_vec(du2), mass_matrix, _vec(du1))
                 @.. du2 -= du
             end
             EEst = norm(du2) / norm(integrator.opts.abstol .+ integrator.opts.reltol .* ks[2])
             integrator.EEst = max(EEst, integrator.EEst)
         end
    end
    cache.linsolve = linres.cache
end

@RosenbrockW6S4OS(:init)
@RosenbrockW6S4OS(:performstep)

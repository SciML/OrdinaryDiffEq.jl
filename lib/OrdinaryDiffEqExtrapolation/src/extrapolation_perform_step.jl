function initialize!(integrator, cache::AitkenNevilleCache)
    integrator.kshortsize = 2
    @unpack k, fsalfirst = cache

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    cache.step_no = 1
    alg = unwrap_alg(integrator, false)
    cache.cur_order = max(alg.init_order, alg.min_order)
end

function perform_step!(integrator, cache::AitkenNevilleCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack k, fsalfirst, T, utilde, atmp, dtpropose, cur_order, A = cache
    @unpack u_tmps, k_tmps = cache

    max_order = min(size(T, 1), cur_order + 1)

    if !isthreaded(alg.threading)
        for i in 1:max_order
            dt_temp = dt / (2^(i - 1))
            # Solve using Euler method
            @muladd @.. broadcast=false u=uprev + dt_temp * fsalfirst
            f(k, u, p, t + dt_temp)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            for j in 2:(2^(i - 1))
                @muladd @.. broadcast=false u=u + dt_temp * k
                f(k, u, p, t + j * dt_temp)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            end
            @.. broadcast=false T[i, 1]=u
        end
    else
        let max_order = max_order, uprev = uprev, dt = dt, fsalfirst = fsalfirst, p = p,
            t = t,
            u_tmps = u_tmps, k_tmps = k_tmps, T = T
            # Balance workload of threads by computing T[1,1] with T[max_order,1] on
            # same thread, T[2,1] with T[max_order-1,1] on same thread. Similarly fill
            # first column of T matrix
            @threaded alg.threading for i in 1:2
                startIndex = (i == 1) ? 1 : max_order
                endIndex = (i == 1) ? max_order - 1 : max_order
                for index in startIndex:endIndex
                    dt_temp = dt / (2^(index - 1))
                    # Solve using Euler method
                    @muladd @.. broadcast=false u_tmps[Threads.threadid()]=uprev +
                                                                           dt_temp *
                                                                           fsalfirst
                    f(k_tmps[Threads.threadid()], u_tmps[Threads.threadid()], p,
                        t + dt_temp)
                    for j in 2:(2^(index - 1))
                        @muladd @.. broadcast=false u_tmps[Threads.threadid()]=u_tmps[Threads.threadid()] +
                                                                               dt_temp *
                                                                               k_tmps[Threads.threadid()]
                        f(k_tmps[Threads.threadid()], u_tmps[Threads.threadid()], p,
                            t + j * dt_temp)
                    end
                    @.. broadcast=false T[index, 1]=u_tmps[Threads.threadid()]
                end
            end
        end
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)^max_order - 1
    end

    # Richardson extrapolation
    tmp = 1
    for j in 2:max_order
        tmp *= 2
        for i in j:max_order
            @.. broadcast=false T[i, j]=(tmp * T[i, j - 1] - T[i - 1, j - 1]) / (tmp - 1)
        end
    end

    if integrator.opts.adaptive
        minimum_work = Inf
        if isone(cache.step_no)
            range_start = 2
        else
            range_start = max(2, cur_order - 1)
        end

        for i in range_start:max_order
            A = 2^(i - 1)
            @.. broadcast=false utilde=T[i, i] - T[i, i - 1]
            atmp = calculate_residuals(utilde, uprev, T[i, i], integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            EEst = integrator.opts.internalnorm(atmp, t)

            beta1 = integrator.opts.controller.beta1
            e = integrator.EEst
            qold = integrator.qold

            integrator.opts.controller.beta1 = 1 / (i + 1)
            integrator.EEst = EEst
            dtpropose = step_accept_controller!(integrator, alg,
                stepsize_controller!(integrator, alg))
            integrator.EEst = e
            integrator.opts.controller.beta1 = beta1
            integrator.qold = qold

            work = A / dtpropose

            if work < minimum_work
                integrator.opts.controller.beta1 = 1 / (i + 1)
                cache.dtpropose = dtpropose
                cache.cur_order = i
                minimum_work = work
                integrator.EEst = EEst
            end
        end
    end

    # using extrapolated value of u
    @.. broadcast=false u=T[cache.cur_order, cache.cur_order]
    cache.step_no = cache.step_no + 1
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::AitkenNevilleConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    cache.step_no = 1
    alg = unwrap_alg(integrator, false)
    cache.cur_order = max(alg.init_order, alg.min_order)
end

function perform_step!(integrator, cache::AitkenNevilleConstantCache, repeat_step = false)
    @unpack t, dt, uprev, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack dtpropose, T, cur_order, work, A = cache

    max_order = min(size(T, 1), cur_order + 1)

    if !isthreaded(alg.threading)
        for i in 1:max_order
            dt_temp = dt / (2^(i - 1)) # Romberg sequence

            # Solve using Euler method with dt_temp = dt/n_{i}
            @muladd u = @.. broadcast=false uprev+dt_temp * integrator.fsalfirst
            k = f(u, p, t + dt_temp)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

            for j in 2:(2^(i - 1))
                @muladd u = @.. broadcast=false u+dt_temp * k
                k = f(u, p, t + j * dt_temp)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            end
            T[i, 1] = u
        end
    else
        let max_order = max_order, dt = dt, uprev = uprev, integrator = integrator, p = p,
            t = t, T = T
            # Balance workload of threads by computing T[1,1] with T[max_order,1] on
            # same thread, T[2,1] with T[max_order-1,1] on same thread. Similarly fill
            # first column of T matrix
            @threaded alg.threading for i in 1:2
                startIndex = (i == 1) ? 1 : max_order
                endIndex = (i == 1) ? max_order - 1 : max_order

                for index in startIndex:endIndex
                    dt_temp = dt / 2^(index - 1)
                    @muladd u = @.. broadcast=false uprev+dt_temp * integrator.fsalfirst
                    k_temp = f(u, p, t + dt_temp)
                    for j in 2:(2^(index - 1))
                        @muladd u = @.. broadcast=false u+dt_temp * k_temp
                        k_temp = f(u, p, t + j * dt_temp)
                    end
                    T[index, 1] = u
                end
            end
        end

        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)^max_order - 1
    end

    # Richardson extrapolation
    tmp = 1
    for j in 2:max_order
        tmp *= 2
        for i in j:max_order
            T[i, j] = (tmp * T[i, j - 1] - T[i - 1, j - 1]) / (tmp - 1)
        end
    end

    if integrator.opts.adaptive
        minimum_work = Inf
        if isone(cache.step_no)
            range_start = 2
        else
            range_start = max(2, cur_order - 1)
        end

        for i in range_start:max_order
            A = 2^(i - 1)
            utilde = T[i, i] - T[i, i - 1]
            atmp = calculate_residuals(utilde, uprev, T[i, i], integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            EEst = integrator.opts.internalnorm(atmp, t)

            beta1 = integrator.opts.controller.beta1
            e = integrator.EEst
            qold = integrator.qold

            integrator.opts.controller.beta1 = 1 / (i + 1)
            integrator.EEst = EEst
            dtpropose = step_accept_controller!(integrator, alg,
                stepsize_controller!(integrator, alg))
            integrator.EEst = e
            integrator.opts.controller.beta1 = beta1
            integrator.qold = qold

            work = A / dtpropose

            if work < minimum_work
                integrator.opts.controller.beta1 = 1 / (i + 1)
                cache.dtpropose = dtpropose
                cache.cur_order = i
                minimum_work = work
                integrator.EEst = EEst
            end
        end
    end

    cache.step_no = cache.step_no + 1

    # Use extrapolated value of u
    integrator.u = T[cache.cur_order, cache.cur_order]

    k = f(integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = k
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitEulerExtrapolationCache)
    integrator.kshortsize = 2
    integrator.f(integrator.fsalfirst, integrator.u, integrator.p, integrator.t)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    cache.step_no = 1
    #alg = unwrap_alg(integrator, true)
    #cache.cur_order = max(alg.init_order, alg.min_order)
end

function perform_step!(integrator, cache::ImplicitEulerExtrapolationCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack T, utilde, atmp, dtpropose, n_curr, A, stage_number, diff1, diff2 = cache
    @unpack J, W, uf, tf, jac_config = cache
    @unpack u_tmps, k_tmps, linsolve_tmps, u_tmps2 = cache

    @unpack sequence = cache

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    if !isthreaded(alg.threading)
        calc_J!(J, integrator, cache) # Store the calculated jac as it won't change in internal discretisation
        for index in 1:(n_curr + 1)
            dt_temp = dt / sequence[index]
            jacobian2W!(W[1], integrator.f.mass_matrix, dt_temp, J)
            integrator.stats.nw += 1
            @.. broadcast=false k_tmps[1]=integrator.fsalfirst
            @.. broadcast=false u_tmps[1]=uprev

            for j in 1:sequence[index]
                @.. broadcast=false linsolve_tmps[1]=k_tmps[1]

                linsolve = cache.linsolve[1]

                if !repeat_step && j == 1
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k_tmps[1]))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k_tmps[1]))
                end

                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false u_tmps2[1]=u_tmps[1]
                @.. broadcast=false u_tmps[1]=u_tmps[1] - k_tmps[1]
                if index <= 2 && j >= 2
                    # Deuflhard Stability check for initial two sequences
                    @.. broadcast=false diff2[1]=u_tmps[1] - u_tmps2[1]
                    @.. broadcast=false diff2[1]=0.5 * (diff2[1] - diff1[1])
                    if integrator.opts.internalnorm(diff1[1], t) <
                       integrator.opts.internalnorm(diff2[1], t)
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                @.. broadcast=false diff1[1]=u_tmps[1] - u_tmps2[1]

                f(k_tmps[1], u_tmps[1], p, t + j * dt_temp)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            end

            @.. broadcast=false T[index, 1]=u_tmps[1]
        end
    else
        calc_J!(J, integrator, cache) # Store the calculated jac as it won't change in internal discretisation
        let n_curr = n_curr, uprev = uprev, dt = dt, p = p, t = t, T = T, W = W,
            integrator = integrator, cache = cache, repeat_step = repeat_step,
            k_tmps = k_tmps, u_tmps = u_tmps, u_tmps2 = u_tmps2, diff1 = diff1,
            diff2 = diff2

            @threaded alg.threading for i in 1:2
                startIndex = (i == 1) ? 1 : n_curr + 1
                endIndex = (i == 1) ? n_curr : n_curr + 1
                for index in startIndex:endIndex
                    dt_temp = dt / sequence[index]
                    jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix, dt_temp, J)
                    @.. broadcast=false k_tmps[Threads.threadid()]=integrator.fsalfirst
                    @.. broadcast=false u_tmps[Threads.threadid()]=uprev
                    for j in 1:sequence[index]
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()]

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step && j == 1
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end

                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false u_tmps2[Threads.threadid()]=u_tmps[Threads.threadid()]
                        @.. broadcast=false u_tmps[Threads.threadid()]=u_tmps[Threads.threadid()] -
                                                                       k_tmps[Threads.threadid()]
                        if index <= 2 && j >= 2
                            # Deuflhard Stability check for initial two sequences
                            @.. broadcast=false diff2[Threads.threadid()]=u_tmps[Threads.threadid()] -
                                                                          u_tmps2[Threads.threadid()]
                            @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                          (diff2[Threads.threadid()] -
                                                                           diff1[Threads.threadid()])
                            if integrator.opts.internalnorm(diff1[Threads.threadid()], t) <
                               integrator.opts.internalnorm(diff2[Threads.threadid()], t)
                                # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                integrator.force_stepfail = true
                                return
                            end
                        end
                        @.. broadcast=false diff1[Threads.threadid()]=u_tmps[Threads.threadid()] -
                                                                      u_tmps2[Threads.threadid()]
                        f(k_tmps[Threads.threadid()], u_tmps[Threads.threadid()], p,
                            t + j * dt_temp)
                    end

                    @.. broadcast=false T[index, 1]=u_tmps[Threads.threadid()]
                end
                integrator.force_stepfail ? break : continue
            end
        end

        nevals = sum(sequence[1:(n_curr + 1)]) - 1
        integrator.stats.nw += n_curr + 1
        integrator.stats.nf += nevals
        integrator.stats.nsolve += nevals
    end

    if integrator.force_stepfail
        return
    end

    # Polynomial extrapolation
    for j in 2:(n_curr + 1)
        for i in j:(n_curr + 1)
            @.. broadcast=false T[i, j]=((sequence[i] / sequence[i - j + 1]) * T[i, j - 1] -
                                         T[i - 1, j - 1]) /
                                        ((sequence[i] / sequence[i - j + 1]) - 1)
        end
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min
            @.. broadcast=false integrator.u=T[i + 1, i + 1]
            @.. broadcast=false cache.utilde=T[i + 1, i]

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(sequence[(n_curr + 2):(win_max + 1)] .//
                                                sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                dt_temp = dt / sequence[n_curr + 1]
                jacobian2W!(W[1], integrator.f.mass_matrix, dt_temp, J)
                integrator.stats.nw += 1
                @.. broadcast=false k_tmps[1]=integrator.fsalfirst
                @.. broadcast=false u_tmps[1]=uprev

                for j in 1:sequence[n_curr + 1]
                    @.. broadcast=false linsolve_tmps[1]=k_tmps[1]

                    linsolve = cache.linsolve[1]

                    if !repeat_step && j == 1
                        linres = dolinsolve(integrator, linsolve; A = W[1],
                            b = _vec(linsolve_tmps[1]),
                            linu = _vec(k_tmps[1]))
                    else
                        linres = dolinsolve(integrator, linsolve; A = nothing,
                            b = _vec(linsolve_tmps[1]),
                            linu = _vec(k_tmps[1]))
                    end

                    cache.linsolve[1] = linres.cache

                    integrator.stats.nsolve += 1
                    @.. broadcast=false u_tmps[1]=u_tmps[1] - k_tmps[1]
                    f(k_tmps[1], u_tmps[1], p, t + j * dt_temp)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                end

                @.. broadcast=false T[n_curr + 1, 1]=u_tmps[1]

                for j in 2:(n_curr + 1)
                    for i in j:(n_curr + 1)
                        @.. broadcast=false T[i, j]=((sequence[i] / sequence[i - j + 1]) *
                                                     T[i, j - 1] - T[i - 1, j - 1]) /
                                                    ((sequence[i] / sequence[i - j + 1]) -
                                                     1)
                    end
                end

                @.. broadcast=false integrator.u=T[n_curr + 1, n_curr + 1]
                @.. broadcast=false cache.utilde=T[n_curr + 1, n_curr]

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        @.. broadcast=false integrator.u=T[n_curr + 1, n_curr + 1]
    end

    cache.step_no = cache.step_no + 1
    f(integrator.fsallast, integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ImplicitEulerExtrapolationConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::ImplicitEulerExtrapolationConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack dtpropose, T, n_curr, work, A, tf, uf = cache
    @unpack sequence, stage_number = cache

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    J = calc_J(integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for index in 1:(n_curr + 1)
            dt_temp = dt / sequence[index]
            W = dt_temp * J - integrator.f.mass_matrix
            integrator.stats.nw += 1
            k_copy = integrator.fsalfirst
            u_tmp = uprev
            diff1 = zero(u_tmp)
            for j in 1:sequence[index]
                k = _reshape(W \ -_vec(dt_temp * k_copy), axes(uprev))
                integrator.stats.nsolve += 1
                u_tmp2 = u_tmp
                u_tmp = u_tmp + k
                if index <= 2 && j >= 2
                    # Deuflhard Stability check for initial two sequences
                    diff2 = u_tmp - u_tmp2
                    diff2 = 0.5 * (diff2 - diff1)
                    if integrator.opts.internalnorm(diff1, t) <
                       integrator.opts.internalnorm(diff2, t)
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                diff1 = u_tmp - u_tmp2
                k_copy = f(u_tmp, p, t + j * dt_temp)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            end

            T[index, 1] = u_tmp
        end
    else
        J = calc_J(integrator, cache) # Store the calculated jac as it won't change in internal discretisation
        let n_curr = n_curr, dt = dt, integrator = integrator, cache = cache,
            repeat_step = repeat_step,
            uprev = uprev, T = T

            @threaded alg.threading for i in 1:2
                startIndex = (i == 1) ? 1 : n_curr + 1
                endIndex = (i == 1) ? n_curr : n_curr + 1
                for index in startIndex:endIndex
                    dt_temp = dt / sequence[index]
                    W = dt_temp * J - integrator.f.mass_matrix
                    k_copy = integrator.fsalfirst
                    u_tmp = uprev
                    diff1 = zero(u_tmp)
                    for j in 1:sequence[index]
                        k = _reshape(W \ -_vec(dt_temp * k_copy), axes(uprev))
                        u_tmp2 = u_tmp
                        u_tmp = u_tmp + k
                        if index <= 2 && j >= 2
                            # Deuflhard Stability check for initial two sequences
                            diff2 = u_tmp - u_tmp2
                            diff2 = 0.5 * (diff2 - diff1)
                            if integrator.opts.internalnorm(diff1, t) <
                               integrator.opts.internalnorm(diff2, t)
                                # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                integrator.force_stepfail = true
                                return
                            end
                        end
                        diff1 = u_tmp - u_tmp2
                        k_copy = f(u_tmp, p, t + j * dt_temp)
                    end
                    T[index, 1] = u_tmp
                end
                integrator.force_stepfail ? break : continue
            end
        end

        if integrator.force_stepfail
            return
        end

        nevals = sum(sequence[1:(n_curr + 1)]) - 1
        integrator.stats.nw += n_curr + 1
        integrator.stats.nf += nevals
        integrator.stats.nsolve += nevals
    end

    # Richardson extrapolation
    tmp = 1
    for j in 2:(n_curr + 1)
        tmp *= 2
        for i in j:(n_curr + 1)
            T[i, j] = ((sequence[i] / sequence[i - j + 1]) * T[i, j - 1] -
                       T[i - 1, j - 1]) /
                      ((sequence[i] / sequence[i - j + 1]) - 1)
        end
    end

    integrator.dt = dt

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min
            u = T[i + 1, i + 1]
            utilde = T[i + 1, i]
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(sequence[(n_curr + 2):(win_max + 1)] .//
                                                sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Always compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                dt_temp = dt / sequence[n_curr + 1]
                W = dt_temp * J - integrator.f.mass_matrix
                integrator.stats.nw += 1
                k_copy = integrator.fsalfirst
                u_tmp = uprev

                for j in 1:sequence[n_curr + 1]
                    k = _reshape(W \ -_vec(dt_temp * k_copy), axes(uprev))
                    integrator.stats.nsolve += 1
                    u_tmp = u_tmp + k
                    k_copy = f(u_tmp, p, t + j * dt_temp)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                end

                T[n_curr + 1, 1] = u_tmp

                #Extrapolate to new order
                for j in 2:(n_curr + 1)
                    for i in j:(n_curr + 1)
                        T[i, j] = ((sequence[i] / sequence[i - j + 1]) * T[i, j - 1] -
                                   T[i - 1, j - 1]) /
                                  ((sequence[i] / sequence[i - j + 1]) - 1)
                    end
                end
                # Update u, integrator.EEst and cache.Q
                u = T[n_curr + 1, n_curr + 1]
                utilde = T[n_curr + 1, n_curr]
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        integrator.u = T[n_curr + 1, n_curr + 1]
    end

    # Use extrapolated value of u
    integrator.u = T[n_curr + 1, n_curr + 1]
    k_temp = f(integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = k_temp
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ExtrapolationMidpointDeuflhardCache)
    # cf. initialize! of MidpointCache
    @unpack k, fsalfirst = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ExtrapolationMidpointDeuflhardCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack n_curr, u_temp1, u_temp2, utilde, res, T, fsalfirst, k = cache
    @unpack u_temp3, u_temp4, k_tmps = cache

    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack stage_number = cache
    @unpack sequence_factor = alg

    fill!(cache.Q, zero(eltype(cache.Q)))
    tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Used by the convergence monitor

    if integrator.opts.adaptive
        # Set up the order window
        win_min = max(alg.min_order, n_curr - 1)
        win_max = min(alg.max_order, n_curr + 1)

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            @.. broadcast=false u_temp2=uprev
            @.. broadcast=false u_temp1=u_temp2 + dt_int * fsalfirst # Euler starting step
            for j in 2:j_int
                f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast=false T[i + 1]=u_temp2 + 2 * dt_int * k # Explicit Midpoint rule
                @.. broadcast=false u_temp2=u_temp1
                @.. broadcast=false u_temp1=T[i + 1]
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr
                    for index in startIndex:endIndex
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        dt_int_temp *
                                                                        fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false T[index + 1]=u_temp4[Threads.threadid()] +
                                                             2 * dt_int_temp *
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                        end
                    end
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = (i, n_curr - i)
                    for index in indices
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        dt_int_temp *
                                                                        fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()], u_temp3[Threads.threadid()], p,
                                t + (j - 1) * dt_int_temp)
                            @.. broadcast=false T[index + 1]=u_temp4[Threads.threadid()] +
                                                             2 * dt_int_temp *
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                        end
                        if indices[2] <= indices[1]
                            break
                        end
                    end
                end
            end
        end
        nevals = cache.stage_number[n_curr - alg.min_order + 1] - 1
        integrator.stats.nf += nevals
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (alg.min_order):n_curr

            #integrator.u .= extrapolation_scalars[i+1] * sum( broadcast(*, cache.T[1:(i+1)], extrapolation_weights[1:(i+1), (i+1)]) ) # Approximation of extrapolation order i
            #cache.utilde .= extrapolation_scalars_2[i] * sum( broadcast(*, cache.T[2:(i+1)], extrapolation_weights_2[1:i, i]) ) # and its internal counterpart

            u_temp1 .= false
            u_temp2 .= false
            for j in 1:(i + 1)
                @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (i + 1)]
            end
            for j in 2:(i + 1)
                @.. broadcast=false u_temp2+=cache.T[j] * extrapolation_weights_2[j - 1, i]
            end
            @.. broadcast=false integrator.u=extrapolation_scalars[i + 1] * u_temp1
            @.. broadcast=false cache.utilde=extrapolation_scalars_2[i] * u_temp2

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif integrator.EEst <=
                   tol^(stage_number[n_curr - alg.min_order + 1] /
                        stage_number[win_max - alg.min_order + 1] - 1)
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update cache.T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                @.. broadcast=false u_temp2=uprev
                @.. broadcast=false u_temp1=u_temp2 + dt_int * fsalfirst # Euler starting step
                for j in 2:j_int
                    f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    @.. broadcast=false T[n_curr + 1]=u_temp2 + 2 * dt_int * k
                    @.. broadcast=false u_temp2=u_temp1
                    @.. broadcast=false u_temp1=T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
                #cache.utilde .= extrapolation_scalars_2[n_curr] * sum( broadcast(*, cache.T[2:(n_curr+1)], extrapolation_weights_2[1:n_curr, n_curr]) ) # and its internal counterpart

                u_temp1 .= false
                u_temp2 .= false
                for j in 1:(n_curr + 1)
                    @.. broadcast=false u_temp1+=cache.T[j] *
                                                 extrapolation_weights[j, (n_curr + 1)]
                end
                for j in 2:(n_curr + 1)
                    @.. broadcast=false u_temp2+=cache.T[j] *
                                                 extrapolation_weights_2[j - 1, n_curr]
                end
                @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
                @.. broadcast=false cache.utilde=extrapolation_scalars_2[n_curr] * u_temp2

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else

        #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
        u_temp1 .= false
        for j in 1:(n_curr + 1)
            @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (n_curr + 1)]
        end
        @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
    end

    f(cache.k, integrator.u, p, t + dt) # Update FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ExtrapolationMidpointDeuflhardConstantCache)
    # cf. initialize! of MidpointConstantCache
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::ExtrapolationMidpointDeuflhardConstantCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack n_curr = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack stage_number = cache
    @unpack sequence_factor = alg

    # Create auxiliary variables
    u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
    u, utilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
    tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Used by the convergence monitor
    T = fill(zero(uprev), alg.max_order + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
    fill!(cache.Q, zero(eltype(cache.Q)))

    # Start computation
    if integrator.opts.adaptive
        # Set up the order window
        win_min = max(alg.min_order, n_curr - 1)
        win_max = min(alg.max_order, n_curr + 1)

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    # Compute the internal discretisations
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            u_temp2 = uprev
            u_temp1 = u_temp2 + dt_int * integrator.fsalfirst # Euler starting step
            for j in 2:j_int
                T[i + 1] = u_temp2 + 2 * dt_int * f(u_temp1, p, t + (j - 1) * dt_int) # Explicit Midpoint rule
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                u_temp2 = u_temp1
                u_temp1 = T[i + 1]
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt,
                integrator = integrator, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr
                    for index in startIndex:endIndex
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        u_temp4 = uprev
                        u_temp3 = u_temp4 + dt_int_temp * integrator.fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            T[index + 1] = u_temp4 +
                                           2 * dt_int_temp *
                                           f(u_temp3, p, t + (j - 1) * dt_int_temp) # Explicit Midpoint rule
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                        end
                    end
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, dt = dt,
                uprev = uprev,
                p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = (i, n_curr - i)
                    for index in indices
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        u_temp4 = uprev
                        u_temp3 = u_temp4 + dt_int_temp * integrator.fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            T[index + 1] = u_temp4 +
                                           2 * dt_int_temp *
                                           f(u_temp3, p, t + (j - 1) * dt_int_temp) # Explicit Midpoint rule
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                        end
                        if indices[2] <= indices[1]
                            break
                        end
                    end
                end
            end
        end
        nevals = cache.stage_number[n_curr - alg.min_order + 1] - 1
        integrator.stats.nf += nevals
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (alg.min_order):n_curr
            u = eltype(uprev).(extrapolation_scalars[i + 1]) *
                sum(broadcast(*, T[1:(i + 1)],
                eltype(uprev).(extrapolation_weights[1:(i + 1), (i + 1)]))) # Approximation of extrapolation order i
            utilde = eltype(uprev).(extrapolation_scalars_2[i]) *
                     sum(broadcast(*, T[2:(i + 1)],
                eltype(uprev).(extrapolation_weights_2[1:i, i]))) # and its internal counterpart
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif integrator.EEst <=
                   tol^(stage_number[n_curr - alg.min_order + 1] /
                        stage_number[win_max - alg.min_order + 1] - 1)
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                u_temp2 = uprev
                u_temp1 = u_temp2 + dt_int * integrator.fsalfirst # Euler starting step
                for j in 2:j_int
                    T[n_curr + 1] = u_temp2 +
                                    2 * dt_int * f(u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    u_temp2 = u_temp1
                    u_temp1 = T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
                    sum(broadcast(*, T[1:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                        (n_curr + 1)]))) # Approximation of extrapolation order n_curr
                utilde = eltype(uprev).(extrapolation_scalars_2[n_curr]) *
                         sum(broadcast(*, T[2:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights_2[1:n_curr,
                        n_curr]))) # and its internal counterpart
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
            sum(broadcast(*, T[1:(n_curr + 1)],
            eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                (n_curr + 1)]))) # Approximation of extrapolation order n_curr
    end

    # Save the latest approximation and update FSAL
    integrator.u = u
    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitDeuflhardExtrapolationCache)
    # cf. initialize! of MidpointCache
    @unpack k, fsalfirst = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ImplicitDeuflhardExtrapolationCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr, u_temp1, u_temp2, utilde, res, T, fsalfirst, k, diff1, diff2 = cache
    @unpack u_temp3, u_temp4, k_tmps = cache

    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack stage_number = cache

    @unpack J, W, uf, tf, linsolve_tmps, jac_config = cache

    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        win_min = max(alg.min_order, n_curr - 1)
        win_max = min(alg.max_order, n_curr + 1)

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    calc_J!(J, integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = 4 * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
            integrator.stats.nw += 1
            @.. broadcast=false u_temp2=uprev
            @.. broadcast=false linsolve_tmps[1]=fsalfirst

            linsolve = cache.linsolve[1]

            if !repeat_step
                linres = dolinsolve(integrator, linsolve; A = W[1],
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            else
                linres = dolinsolve(integrator, linsolve; A = nothing,
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            end

            cache.linsolve[1] = linres.cache

            integrator.stats.nsolve += 1
            @.. broadcast=false u_temp1=u_temp2 - k # Euler starting step
            @.. broadcast=false diff1[1]=u_temp1 - u_temp2
            for j in 2:j_int
                f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast=false linsolve_tmps[1]=k - (u_temp1 - u_temp2) / dt_int

                linsolve = cache.linsolve[1]

                if !repeat_step && j == 1
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                end
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false T[i + 1]=2 * u_temp1 - u_temp2 - 2 * k # Explicit Midpoint rule
                @.. broadcast=false u_temp2=u_temp1
                @.. broadcast=false u_temp1=T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    @.. broadcast=false diff2[1]=u_temp1 - u_temp2
                    if (integrator.opts.internalnorm(diff1[1], t) <
                        integrator.opts.internalnorm(0.5 * (diff2[1] - diff1[1]), t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int_temp = 4 * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix,
                            dt_int_temp, J)
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=fsalfirst

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end
                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false k_tmps[Threads.threadid()]=-k_tmps[Threads.threadid()]
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        k_tmps[Threads.threadid()] # Euler starting step
                        @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                      u_temp4[Threads.threadid()]
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()] -
                                                                                  (u_temp3[Threads.threadid()] -
                                                                                   u_temp4[Threads.threadid()]) /
                                                                                  dt_int_temp

                            linsolve = cache.linsolve[Threads.threadid()]

                            if !repeat_step && j == 1
                                linres = dolinsolve(integrator, linsolve;
                                    A = W[Threads.threadid()],
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            end
                            cache.linsolve[Threads.threadid()] = linres.cache

                            @.. broadcast=false T[index + 1]=2 *
                                                             u_temp3[Threads.threadid()] -
                                                             u_temp4[Threads.threadid()] -
                                                             2 * k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                              u_temp4[Threads.threadid()]
                                @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                              (diff2[Threads.threadid()] -
                                                                               diff1[Threads.threadid()])
                                if (integrator.opts.internalnorm(diff1[Threads.threadid()],
                                    t) <
                                    integrator.opts.internalnorm(diff2[Threads.threadid()],
                                    t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i) #Use flag to avoid union
                    for index in indices
                        index == -1 && continue
                        j_int_temp = 4 * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix,
                            dt_int_temp, J)
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=fsalfirst

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end
                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false k_tmps[Threads.threadid()]=-k_tmps[Threads.threadid()]
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        k_tmps[Threads.threadid()] # Euler starting step
                        @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                      u_temp4[Threads.threadid()]
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()] -
                                                                                  (u_temp3[Threads.threadid()] -
                                                                                   u_temp4[Threads.threadid()]) /
                                                                                  dt_int_temp

                            linsolve = cache.linsolve[Threads.threadid()]

                            if !repeat_step && j == 1
                                linres = dolinsolve(integrator, linsolve;
                                    A = W[Threads.threadid()],
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            end
                            cache.linsolve[Threads.threadid()] = linres.cache

                            @.. broadcast=false T[index + 1]=2 *
                                                             u_temp3[Threads.threadid()] -
                                                             u_temp4[Threads.threadid()] -
                                                             2 * k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                              u_temp4[Threads.threadid()]
                                @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                              (diff2[Threads.threadid()] -
                                                                               diff1[Threads.threadid()])
                                if (integrator.opts.internalnorm(diff1[Threads.threadid()],
                                    t) <
                                    integrator.opts.internalnorm(diff2[Threads.threadid()],
                                    t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (alg.min_order):n_curr

            #integrator.u .= extrapolation_scalars[i+1] * sum( broadcast(*, cache.T[1:(i+1)], extrapolation_weights[1:(i+1), (i+1)]) ) # Approximation of extrapolation order i
            #cache.utilde .= extrapolation_scalars_2[i] * sum( broadcast(*, cache.T[2:(i+1)], extrapolation_weights_2[1:i, i]) ) # and its internal counterpart

            u_temp1 .= false
            u_temp2 .= false
            for j in 1:(i + 1)
                @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (i + 1)]
            end
            for j in 2:(i + 1)
                @.. broadcast=false u_temp2+=cache.T[j] * extrapolation_weights_2[j - 1, i]
            end
            @.. broadcast=false integrator.u=extrapolation_scalars[i + 1] * u_temp1
            @.. broadcast=false cache.utilde=extrapolation_scalars_2[i] * u_temp2

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        while n_curr <= win_max
            tol = integrator.opts.internalnorm(cache.utilde - integrator.u, t) /
                  integrator.EEst # Used by the convergence monitor
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif integrator.EEst <=
                   tol^(stage_number[n_curr - alg.min_order + 1] /
                        stage_number[win_max - alg.min_order + 1] - 1)
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update cache.T
                j_int = 4 * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
                integrator.stats.nw += 1
                @.. broadcast=false u_temp2=uprev
                @.. broadcast=false linsolve_tmps[1]=fsalfirst

                linsolve = cache.linsolve[1]
                linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmps[1]),
                    linu = _vec(k))
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false u_temp1=u_temp2 - k # Euler starting step
                for j in 2:j_int
                    f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    @.. broadcast=false linsolve_tmps[1]=dt_int * k - (u_temp1 - u_temp2)

                    linsolve = cache.linsolve[1]
                    linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmps[1]),
                        linu = _vec(k))
                    cache.linsolve[1] = linres.cache

                    integrator.stats.nsolve += 1
                    @.. broadcast=false T[n_curr + 1]=2 * u_temp1 - u_temp2 - 2 * k # Explicit Midpoint rule
                    @.. broadcast=false u_temp2=u_temp1
                    @.. broadcast=false u_temp1=T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
                #cache.utilde .= extrapolation_scalars_2[n_curr] * sum( broadcast(*, cache.T[2:(n_curr+1)], extrapolation_weights_2[1:n_curr, n_curr]) ) # and its internal counterpart

                u_temp1 .= false
                u_temp2 .= false
                for j in 1:(n_curr + 1)
                    @.. broadcast=false u_temp1+=cache.T[j] *
                                                 extrapolation_weights[j, (n_curr + 1)]
                end
                for j in 2:(n_curr + 1)
                    @.. broadcast=false u_temp2+=cache.T[j] *
                                                 extrapolation_weights_2[j - 1, n_curr]
                end
                @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
                @.. broadcast=false cache.utilde=extrapolation_scalars_2[n_curr] * u_temp2

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else

        #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
        u_temp1 .= false
        for j in 1:(n_curr + 1)
            @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (n_curr + 1)]
        end
        @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
    end

    f(cache.k, integrator.u, p, t + dt) # Update FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ImplicitDeuflhardExtrapolationConstantCache)
    # cf. initialize! of MidpointConstantCache
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::ImplicitDeuflhardExtrapolationConstantCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack stage_number = cache

    # Create auxiliary variables
    u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
    u, utilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
    T = fill(zero(uprev), alg.max_order + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
    fill!(cache.Q, zero(eltype(cache.Q)))

    # Start computation
    if integrator.opts.adaptive
        # Set up the order window
        win_min = max(alg.min_order, n_curr - 1)
        win_max = min(alg.max_order, n_curr + 1)

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    # Compute the internal discretisations
    J = calc_J(integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = 4 * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            W = dt_int * J - integrator.f.mass_matrix
            integrator.stats.nw += 1
            u_temp2 = uprev
            u_temp1 = u_temp2 +
                      _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
            diff1 = u_temp1 - u_temp2
            for j in 2:j_int
                T[i + 1] = 2 * u_temp1 - u_temp2 +
                           2 * _reshape(
                    W \
                    -_vec(dt_int * f(u_temp1, p, t + (j - 1) * dt_int) -
                          (u_temp1 - u_temp2)),
                    axes(uprev))
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                u_temp2 = u_temp1
                u_temp1 = T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    diff2 = u_temp1 - u_temp2
                    if (integrator.opts.internalnorm(diff1, t) <
                        integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp2 = u_temp2,
                u_temp2 = u_temp2, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int = 4 * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:j_int
                            T[index + 1] = 2 * u_temp3 - u_temp4 +
                                           2 * _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int) -
                                      (u_temp3 - u_temp4)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1[1], t) <
                                    integrator.opts.internalnorm(
                                    0.5 *
                                    (diff2[1] - diff1[1]), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt,
                integrator = integrator, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int = 4 * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:j_int
                            T[index + 1] = 2 * u_temp3 - u_temp4 +
                                           2 * _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int) -
                                      (u_temp3 - u_temp4)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1, t) <
                                    integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (alg.min_order):n_curr
            u = eltype(uprev).(extrapolation_scalars[i + 1]) *
                sum(broadcast(*, T[1:(i + 1)],
                eltype(uprev).(extrapolation_weights[1:(i + 1), (i + 1)]))) # Approximation of extrapolation order i
            utilde = eltype(uprev).(extrapolation_scalars_2[i]) *
                     sum(broadcast(*, T[2:(i + 1)],
                eltype(uprev).(extrapolation_weights_2[1:i, i]))) # and its internal counterpart
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        while n_curr <= win_max
            tol = integrator.opts.internalnorm(utilde - u, t) / integrator.EEst # Used by the convergence monitor
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif integrator.EEst <=
                   tol^(stage_number[n_curr - alg.min_order + 1] /
                        stage_number[win_max - alg.min_order + 1] - 1)
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                j_int = 4 * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                W = dt_int * J - integrator.f.mass_matrix
                integrator.stats.nw += 1
                u_temp2 = uprev
                u_temp1 = u_temp2 +
                          _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
                for j in 2:j_int
                    T[n_curr + 1] = 2 * u_temp1 - u_temp2 +
                                    2 * _reshape(
                        W \
                        -_vec(dt_int *
                              f(u_temp1, p, t + (j - 1) * dt_int) -
                              (u_temp1 - u_temp2)),
                        axes(uprev))
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    u_temp2 = u_temp1
                    u_temp1 = T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
                    sum(broadcast(*, T[1:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                        (n_curr + 1)]))) # Approximation of extrapolation order n_curr
                utilde = eltype(uprev).(extrapolation_scalars_2[n_curr]) *
                         sum(broadcast(*, T[2:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights_2[1:n_curr,
                        n_curr]))) # and its internal counterpart
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
            sum(broadcast(*, T[1:(n_curr + 1)],
            eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                (n_curr + 1)]))) # Approximation of extrapolation order n_curr
    end

    # Save the latest approximation and update FSAL
    integrator.u = u
    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ExtrapolationMidpointHairerWannerCache)
    # cf. initialize! of MidpointCache
    @unpack k, fsalfirst = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ExtrapolationMidpointHairerWannerCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack n_curr, u_temp1, u_temp2, utilde, res, T, fsalfirst, k = cache
    @unpack u_temp3, u_temp4, k_tmps = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack sequence_factor = alg

    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            @.. broadcast=false u_temp2=uprev
            @.. broadcast=false u_temp1=u_temp2 + dt_int * fsalfirst # Euler starting step
            for j in 2:j_int
                f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast=false T[i + 1]=u_temp2 + 2 * dt_int * k # Explicit Midpoint rule
                @.. broadcast=false u_temp2=u_temp1
                @.. broadcast=false u_temp1=T[i + 1]
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        dt_int_temp *
                                                                        fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false T[index + 1]=u_temp4[Threads.threadid()] +
                                                             2 * dt_int_temp *
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                        end
                    end
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] +
                                                                        dt_int_temp *
                                                                        fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false T[index + 1]=u_temp4[Threads.threadid()] +
                                                             2 * dt_int_temp *
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                        end
                    end
                end
            end
        end
        nevals = cache.stage_number[n_curr + 1] - 1
        integrator.stats.nf += nevals
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min

            #integrator.u .= extrapolation_scalars[i+1] * sum( broadcast(*, cache.T[1:(i+1)], extrapolation_weights[1:(i+1), (i+1)]) ) # Approximation of extrapolation order i
            #cache.utilde .= extrapolation_scalars_2[i] * sum( broadcast(*, cache.T[2:(i+1)], extrapolation_weights_2[1:i, i]) ) # and its internal counterpart

            u_temp1 .= false
            u_temp2 .= false
            for j in 1:(i + 1)
                @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (i + 1)]
            end
            for j in 2:(i + 1)
                @.. broadcast=false u_temp2+=cache.T[j] * extrapolation_weights_2[j - 1, i]
            end
            @.. broadcast=false integrator.u=extrapolation_scalars[i + 1] * u_temp1
            @.. broadcast=false cache.utilde=extrapolation_scalars_2[i] * u_temp2

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(subdividing_sequence[(n_curr + 2):(win_max + 1)] .//
                                                subdividing_sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update cache.T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                @.. broadcast=false u_temp2=uprev
                @.. broadcast=false u_temp1=u_temp2 + dt_int * fsalfirst # Euler starting step
                for j in 2:j_int
                    f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    @.. broadcast=false T[n_curr + 1]=u_temp2 + 2 * dt_int * k
                    @.. broadcast=false u_temp2=u_temp1
                    @.. broadcast=false u_temp1=T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
                #cache.utilde .= extrapolation_scalars_2[n_curr] * sum( broadcast(*, cache.T[2:(n_curr+1)], extrapolation_weights_2[1:n_curr, n_curr]) ) # and its internal counterpart

                u_temp1 .= false
                u_temp2 .= false
                for j in 1:(n_curr + 1)
                    @.. broadcast=false u_temp1+=cache.T[j] *
                                                 extrapolation_weights[j, (n_curr + 1)]
                end
                for j in 2:(n_curr + 1)
                    @.. broadcast=false u_temp2+=cache.T[j] *
                                                 extrapolation_weights_2[j - 1, n_curr]
                end
                @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
                @.. broadcast=false cache.utilde=extrapolation_scalars_2[n_curr] * u_temp2

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else

        #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
        u_temp1 .= false
        for j in 1:(n_curr + 1)
            @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (n_curr + 1)]
        end
        @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
    end

    f(cache.k, integrator.u, p, t + dt) # Update FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ExtrapolationMidpointHairerWannerConstantCache)
    # cf. initialize! of MidpointConstantCache
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::ExtrapolationMidpointHairerWannerConstantCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack n_curr = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack sequence_factor = alg

    # Create auxiliary variables
    u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
    u, utilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
    T = fill(zero(uprev), alg.max_order + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            u_temp2 = uprev
            u_temp1 = u_temp2 + dt_int * integrator.fsalfirst # Euler starting step
            for j in 2:j_int
                T[i + 1] = u_temp2 + 2 * dt_int * f(u_temp1, p, t + (j - 1) * dt_int) # Explicit Midpoint rule
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                u_temp2 = u_temp1
                u_temp1 = T[i + 1]
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, dt = dt,
                uprev = uprev,
                integrator = integrator, T = T, p = p, t = t

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr
                    for index in startIndex:endIndex
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        u_temp4 = uprev
                        u_temp3 = u_temp4 + dt_int_temp * integrator.fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            T[index + 1] = u_temp4 +
                                           2 * dt_int_temp *
                                           f(u_temp3, p, t + (j - 1) * dt_int_temp) # Explicit Midpoint rule
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                        end
                    end
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, dt = dt,
                uprev = uprev,
                integrator = integrator, T = T, p = p, t = t

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        u_temp4 = uprev
                        u_temp3 = u_temp4 + dt_int_temp * integrator.fsalfirst # Euler starting step
                        for j in 2:j_int_temp
                            T[index + 1] = u_temp4 +
                                           2 * dt_int_temp *
                                           f(u_temp3, p, t + (j - 1) * dt_int_temp) # Explicit Midpoint rule
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                        end
                    end
                end
            end
        end
        nevals = cache.stage_number[n_curr + 1] - 1
        integrator.stats.nf += nevals
    end
    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min
            u = eltype(uprev).(extrapolation_scalars[i + 1]) *
                sum(broadcast(*, T[1:(i + 1)],
                eltype(uprev).(extrapolation_weights[1:(i + 1), (i + 1)]))) # Approximation of extrapolation order i
            utilde = eltype(uprev).(extrapolation_scalars_2[i]) *
                     sum(broadcast(*, T[2:(i + 1)],
                eltype(uprev).(extrapolation_weights_2[1:i, i]))) # and its internal counterpart
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(subdividing_sequence[(n_curr + 2):(win_max + 1)] .//
                                                subdividing_sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Always compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                u_temp2 = uprev
                u_temp1 = u_temp2 + dt_int * integrator.fsalfirst # Euler starting step
                for j in 2:j_int
                    T[n_curr + 1] = u_temp2 +
                                    2 * dt_int * f(u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    u_temp2 = u_temp1
                    u_temp1 = T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
                    sum(broadcast(*, T[1:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                        (n_curr + 1)]))) # Approximation of extrapolation order n_curr
                utilde = eltype(uprev).(extrapolation_scalars_2[n_curr]) *
                         sum(broadcast(*, T[2:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights_2[1:n_curr,
                        n_curr]))) # and its internal counterpart
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
            sum(broadcast(*, T[1:(n_curr + 1)],
            eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                (n_curr + 1)]))) # Approximation of extrapolation order n_curr
    end

    # Save the latest approximation and update FSAL
    integrator.u = u
    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ImplicitHairerWannerExtrapolationConstantCache)
    # cf. initialize! of MidpointConstantCache
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::ImplicitHairerWannerExtrapolationConstantCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients

    # Create auxiliary variables
    u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
    u, utilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
    T = fill(zero(uprev), alg.max_order + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    J = calc_J(integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = 4 * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            W = dt_int * J - integrator.f.mass_matrix
            integrator.stats.nw += 1
            u_temp2 = uprev
            u_temp1 = u_temp2 +
                      _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
            diff1 = u_temp1 - u_temp2
            for j in 2:(j_int + 1)
                T[i + 1] = 2 * u_temp1 - u_temp2 +
                           2 * _reshape(
                    W \
                    -_vec(dt_int * f(u_temp1, p, t + (j - 1) * dt_int) -
                          (u_temp1 - u_temp2)),
                    axes(uprev))
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                if (j == j_int + 1)
                    T[i + 1] = 0.5(T[i + 1] + u_temp2)
                end
                u_temp2 = u_temp1
                u_temp1 = T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    diff2 = u_temp1 - u_temp2
                    if (integrator.opts.internalnorm(diff1, t) <
                        integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                diff1 = u_temp1 - u_temp2
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp2 = u_temp2,
                u_temp2 = u_temp2, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int = 4 * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:(j_int + 1)
                            T[index + 1] = 2 * u_temp3 - u_temp4 +
                                           2 * _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int) -
                                      (u_temp3 - u_temp4)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            if (j == j_int + 1)
                                T[index + 1] = 0.5(T[index + 1] + u_temp4)
                            end
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1, t) <
                                    integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            diff1 = u_temp3 - u_temp4
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt,
                integrator = integrator, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int = 4 * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:(j_int + 1)
                            T[index + 1] = 2 * u_temp3 - u_temp4 +
                                           2 * _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int) -
                                      (u_temp3 - u_temp4)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            if (j == j_int + 1)
                                T[index + 1] = 0.5(T[index + 1] + u_temp4)
                            end
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1, t) <
                                    integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            diff1 = u_temp3 - u_temp4
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min
            u = eltype(uprev).(extrapolation_scalars[i + 1]) *
                sum(broadcast(*, T[1:(i + 1)],
                eltype(uprev).(extrapolation_weights[1:(i + 1), (i + 1)]))) # Approximation of extrapolation order i
            utilde = eltype(uprev).(extrapolation_scalars_2[i]) *
                     sum(broadcast(*, T[2:(i + 1)],
                eltype(uprev).(extrapolation_weights_2[1:i, i]))) # and its internal counterpart
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(subdividing_sequence[(n_curr + 2):(win_max + 1)] .//
                                                subdividing_sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Always compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                j_int = 4 * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                W = dt_int * J - integrator.f.mass_matrix
                integrator.stats.nw += 1
                u_temp2 = uprev
                u_temp1 = u_temp2 +
                          _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
                for j in 2:(j_int + 1)
                    T[n_curr + 1] = 2 * u_temp1 - u_temp2 +
                                    2 * _reshape(
                        W \
                        -_vec(dt_int *
                              f(u_temp1, p, t + (j - 1) * dt_int) -
                              (u_temp1 - u_temp2)),
                        axes(uprev))
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    if (j == j_int + 1)
                        T[n_curr + 1] = 0.5(T[n_curr + 1] + u_temp2)
                    end
                    u_temp2 = u_temp1
                    u_temp1 = T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
                    sum(broadcast(*, T[1:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                        (n_curr + 1)]))) # Approximation of extrapolation order n_curr
                utilde = eltype(uprev).(extrapolation_scalars_2[n_curr]) *
                         sum(broadcast(*, T[2:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights_2[1:n_curr,
                        n_curr]))) # and its internal counterpart
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
            sum(broadcast(*, T[1:(n_curr + 1)],
            eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                (n_curr + 1)]))) # Approximation of extrapolation order n_curr
    end

    # Save the latest approximation and update FSAL
    integrator.u = u
    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitHairerWannerExtrapolationCache)
    # cf. initialize! of MidpointCache
    @unpack k, fsalfirst = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ImplicitHairerWannerExtrapolationCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr, u_temp1, u_temp2, utilde, res, T, fsalfirst, k, diff1, diff2 = cache
    @unpack u_temp3, u_temp4, k_tmps = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients

    @unpack J, W, uf, tf, linsolve_tmps, jac_config = cache

    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    calc_J!(J, integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = 4 * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
            integrator.stats.nw += 1
            @.. broadcast=false u_temp2=uprev
            @.. broadcast=false linsolve_tmps[1]=fsalfirst

            linsolve = cache.linsolve[1]
            if !repeat_step
                linres = dolinsolve(integrator, linsolve; A = W[1],
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            else
                linres = dolinsolve(integrator, linsolve; A = nothing,
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            end
            cache.linsolve[1] = linres.cache

            integrator.stats.nsolve += 1
            @.. broadcast=false u_temp1=u_temp2 - k # Euler starting step
            @.. broadcast=false diff1[1]=u_temp1 - u_temp2
            for j in 2:(j_int + 1)
                f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast=false linsolve_tmps[1]=k - (u_temp1 - u_temp2) / dt_int

                linsolve = cache.linsolve[1]

                if !repeat_step && j == 1
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                end
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false T[i + 1]=2 * u_temp1 - u_temp2 - 2 * k # Explicit Midpoint rule
                if (j == j_int + 1)
                    @.. broadcast=false T[i + 1]=0.5(T[i + 1] + u_temp2)
                end
                @.. broadcast=false u_temp2=u_temp1
                @.. broadcast=false u_temp1=T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    @.. broadcast=false diff2[1]=u_temp1 - u_temp2
                    @.. broadcast=false diff2[1]=0.5 * (diff2[1] - diff1[1])
                    if (integrator.opts.internalnorm(diff1[1], t) <
                        integrator.opts.internalnorm(diff2[1], t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                @.. broadcast=false diff1[1]=u_temp1 - u_temp2
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int_temp = 4 * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix,
                            dt_int_temp, J)
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=fsalfirst

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end
                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] -
                                                                        k_tmps[Threads.threadid()] # Euler starting step
                        @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                      u_temp4[Threads.threadid()]
                        for j in 2:(j_int_temp + 1)
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()] -
                                                                                  (u_temp3[Threads.threadid()] -
                                                                                   u_temp4[Threads.threadid()]) /
                                                                                  dt_int_temp

                            linsolve = cache.linsolve[Threads.threadid()]
                            if !repeat_step && j == 1
                                linres = dolinsolve(integrator, linsolve;
                                    A = W[Threads.threadid()],
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            end
                            cache.linsolve[Threads.threadid()] = linres.cache

                            @.. broadcast=false T[index + 1]=2 *
                                                             u_temp3[Threads.threadid()] -
                                                             u_temp4[Threads.threadid()] -
                                                             2 * k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            if (j == j_int_temp + 1)
                                @.. broadcast=false T[index + 1]=0.5(T[index + 1] +
                                                                     u_temp4[Threads.threadid()])
                            end
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                              u_temp4[Threads.threadid()]
                                @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                              (diff2[Threads.threadid()] -
                                                                               diff1[Threads.threadid()])
                                if (integrator.opts.internalnorm(diff1[Threads.threadid()],
                                    t) <
                                    integrator.opts.internalnorm(diff2[Threads.threadid()],
                                    t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                          u_temp4[Threads.threadid()]
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    tid = Threads.threadid()
                    linsolvetmp = linsolve_tmps[tid]
                    ktmp = k_tmps[tid]
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i) #Use flag to avoid type union/branch
                    for index in indices
                        index == -1 && continue
                        j_int_temp = 4 * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[tid], integrator.f.mass_matrix, dt_int_temp, J)
                        @.. broadcast=false u_temp4[tid]=uprev
                        @.. broadcast=false linsolvetmp=fsalfirst

                        linsolve = cache.linsolve[tid]
                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve; A = W[tid],
                                b = _vec(linsolvetmp), linu = _vec(ktmp))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolvetmp), linu = _vec(ktmp))
                        end
                        cache.linsolve[tid] = linres.cache

                        @.. broadcast=false u_temp3[tid]=u_temp4[tid] - ktmp # Euler starting step
                        @.. broadcast=false diff1[tid]=u_temp3[tid] - u_temp4[tid]
                        for j in 2:(j_int_temp + 1)
                            f(ktmp, cache.u_temp3[tid], p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolvetmp=ktmp -
                                                            (u_temp3[tid] - u_temp4[tid]) /
                                                            dt_int_temp

                            linsolve = cache.linsolve[tid]

                            if (!repeat_step && j == 1)
                                linres = dolinsolve(integrator, linsolve; A = W[tid],
                                    b = _vec(linsolvetmp),
                                    linu = _vec(ktmp))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolvetmp),
                                    linu = _vec(ktmp))
                            end
                            cache.linsolve[tid] = linres.cache

                            @.. broadcast=false T[index + 1]=2 * u_temp3[tid] -
                                                             u_temp4[tid] - 2 * ktmp # Explicit Midpoint rule
                            if (j == j_int_temp + 1)
                                @.. broadcast=false T[index + 1]=0.5(T[index + 1] +
                                                                     u_temp4[tid])
                            end
                            @.. broadcast=false u_temp4[tid]=u_temp3[tid]
                            @.. broadcast=false u_temp3[tid]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[tid]=u_temp3[tid] - u_temp4[tid]
                                @.. broadcast=false diff2[tid]=0.5 *
                                                               (diff2[tid] - diff1[tid])
                                if (integrator.opts.internalnorm(diff1[tid], t) <
                                    integrator.opts.internalnorm(diff2[tid], t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            @.. broadcast=false diff1[tid]=u_temp3[tid] - u_temp4[tid]
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min

            #integrator.u .= extrapolation_scalars[i+1] * sum( broadcast(*, cache.T[1:(i+1)], extrapolation_weights[1:(i+1), (i+1)]) ) # Approximation of extrapolation order i
            #cache.utilde .= extrapolation_scalars_2[i] * sum( broadcast(*, cache.T[2:(i+1)], extrapolation_weights_2[1:i, i]) ) # and its internal counterpart

            u_temp1 .= false
            u_temp2 .= false
            for j in 1:(i + 1)
                @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (i + 1)]
            end
            for j in 2:(i + 1)
                @.. broadcast=false u_temp2+=cache.T[j] * extrapolation_weights_2[j - 1, i]
            end
            @.. broadcast=false integrator.u=extrapolation_scalars[i + 1] * u_temp1
            @.. broadcast=false cache.utilde=extrapolation_scalars_2[i] * u_temp2

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            EEst1 = one(integrator.EEst)
            for i in (n_curr + 2):(win_max + 1)
                EEst1 *= subdividing_sequence[i] / subdividing_sequence[1]
            end
            EEst1 *= EEst1
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) || integrator.EEst <= EEst1
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update cache.T
                j_int = 4 * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
                integrator.stats.nw += 1
                @.. broadcast=false u_temp2=uprev
                @.. broadcast=false linsolve_tmps[1]=fsalfirst

                linsolve = cache.linsolve[1]

                if !repeat_step
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                end
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false u_temp1=u_temp2 - k # Euler starting step
                for j in 2:(j_int + 1)
                    f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    @.. broadcast=false linsolve_tmps[1]=k - (u_temp1 - u_temp2) / dt_int

                    linsolve = cache.linsolve[1]

                    if !repeat_step && j == 1
                        linres = dolinsolve(integrator, linsolve; A = W[1],
                            b = _vec(linsolve_tmps[1]), linu = _vec(k))
                    else
                        linres = dolinsolve(integrator, linsolve; A = nothing,
                            b = _vec(linsolve_tmps[1]), linu = _vec(k))
                    end
                    cache.linsolve[1] = linres.cache

                    integrator.stats.nsolve += 1
                    @.. broadcast=false T[n_curr + 1]=2 * u_temp1 - u_temp2 - 2 * k # Explicit Midpoint rule
                    if (j == j_int + 1)
                        @.. broadcast=false T[n_curr + 1]=0.5(T[n_curr + 1] + u_temp2)
                    end
                    @.. broadcast=false u_temp2=u_temp1
                    @.. broadcast=false u_temp1=T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
                #cache.utilde .= extrapolation_scalars_2[n_curr] * sum( broadcast(*, cache.T[2:(n_curr+1)], extrapolation_weights_2[1:n_curr, n_curr]) ) # and its internal counterpart

                u_temp1 .= false
                u_temp2 .= false
                for j in 1:(n_curr + 1)
                    @.. broadcast=false u_temp1+=cache.T[j] *
                                                 extrapolation_weights[j, (n_curr + 1)]
                end
                for j in 2:(n_curr + 1)
                    @.. broadcast=false u_temp2+=cache.T[j] *
                                                 extrapolation_weights_2[j - 1, n_curr]
                end
                @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
                @.. broadcast=false cache.utilde=extrapolation_scalars_2[n_curr] * u_temp2

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else

        #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
        u_temp1 .= false
        for j in 1:(n_curr + 1)
            @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (n_curr + 1)]
        end
        @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
    end

    f(cache.k, integrator.u, p, t + dt) # Update FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::ImplicitEulerBarycentricExtrapolationConstantCache)
    # cf. initialize! of MidpointConstantCache
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,
        cache::ImplicitEulerBarycentricExtrapolationConstantCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack sequence_factor = alg

    # Create auxiliary variables
    u_temp1, u_temp2 = copy(uprev), copy(uprev) # Auxiliary variables for computing the internal discretisations
    u, utilde = copy(uprev), copy(uprev) # Storage for the latest approximation and its internal counterpart
    T = fill(zero(uprev), alg.max_order + 1) # Storage for the internal discretisations obtained by the explicit midpoint rule
    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    J = calc_J(integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            W = dt_int * J - integrator.f.mass_matrix
            integrator.stats.nw += 1
            u_temp2 = uprev
            u_temp1 = u_temp2 +
                      _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
            diff1 = u_temp1 - u_temp2
            for j in 2:(j_int + 1)
                T[i + 1] = u_temp1 +
                           _reshape(
                    W \ -_vec(dt_int * f(u_temp1, p, t + (j - 1) * dt_int)),
                    axes(uprev))
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                if (j == j_int + 1)
                    T[i + 1] = 0.25(T[i + 1] + 2 * u_temp1 + u_temp2)
                end
                u_temp2 = u_temp1
                u_temp1 = T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    diff2 = u_temp1 - u_temp2
                    if (integrator.opts.internalnorm(diff1, t) <
                        integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                diff1 = u_temp1 - u_temp2
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp2 = u_temp2,
                u_temp2 = u_temp2, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int = sequence_factor * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:(j_int + 1)
                            T[index + 1] = u_temp3 + _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            if (j == j_int + 1)
                                T[index + 1] = 0.25(T[index + 1] + 2 * u_temp3 + u_temp4)
                            end
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1, t) <
                                    integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            diff1 = u_temp3 - u_temp4
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt,
                integrator = integrator, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int = sequence_factor * subdividing_sequence[index + 1]
                        dt_int = dt / j_int # Stepsize of the ith internal discretisation
                        W = dt_int * J - integrator.f.mass_matrix
                        integrator.stats.nw += 1
                        u_temp4 = uprev
                        u_temp3 = u_temp4 +
                                  _reshape(W \ -_vec(dt_int * integrator.fsalfirst),
                            axes(uprev)) # Euler starting step
                        diff1 = u_temp3 - u_temp4
                        for j in 2:(j_int + 1)
                            T[index + 1] = u_temp3 + _reshape(
                                W \
                                -_vec(dt_int * f(u_temp3, p,
                                    t + (j - 1) * dt_int)),
                                axes(uprev))
                            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                            if (j == j_int + 1)
                                T[index + 1] = 0.25(T[index + 1] + 2 * u_temp3 + u_temp4)
                            end
                            u_temp4 = u_temp3
                            u_temp3 = T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                diff2 = u_temp3 - u_temp4
                                if (integrator.opts.internalnorm(diff1, t) <
                                    integrator.opts.internalnorm(0.5 * (diff2 - diff1), t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            diff1 = u_temp3 - u_temp4
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min
            u = eltype(uprev).(extrapolation_scalars[i + 1]) *
                sum(broadcast(*, T[1:(i + 1)],
                eltype(uprev).(extrapolation_weights[1:(i + 1), (i + 1)]))) # Approximation of extrapolation order i
            utilde = eltype(uprev).(extrapolation_scalars_2[i]) *
                     sum(broadcast(*, T[2:(i + 1)],
                eltype(uprev).(extrapolation_weights_2[1:i, i]))) # and its internal counterpart
            res = calculate_residuals(u, utilde, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) ||
                   integrator.EEst <=
                   typeof(integrator.EEst)(prod(subdividing_sequence[(n_curr + 2):(win_max + 1)] .//
                                                subdividing_sequence[1]^2))
                # Reject current approximation order but pass convergence monitor
                # Always compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                W = dt_int * J - integrator.f.mass_matrix
                integrator.stats.nw += 1
                u_temp2 = uprev
                u_temp1 = u_temp2 +
                          _reshape(W \ -_vec(dt_int * integrator.fsalfirst), axes(uprev)) # Euler starting step
                for j in 2:(j_int + 1)
                    T[n_curr + 1] = u_temp1 + _reshape(
                        W \
                        -_vec(dt_int *
                              f(u_temp1, p, t + (j - 1) * dt_int)),
                        axes(uprev))
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    if (j == j_int + 1)
                        T[n_curr + 1] = 0.25(T[n_curr + 1] + 2 * u_temp1 + u_temp2)
                    end
                    u_temp2 = u_temp1
                    u_temp1 = T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
                    sum(broadcast(*, T[1:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                        (n_curr + 1)]))) # Approximation of extrapolation order n_curr
                utilde = eltype(uprev).(extrapolation_scalars_2[n_curr]) *
                         sum(broadcast(*, T[2:(n_curr + 1)],
                    eltype(uprev).(extrapolation_weights_2[1:n_curr,
                        n_curr]))) # and its internal counterpart
                res = calculate_residuals(u, utilde, integrator.opts.abstol,
                    integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else
        u = eltype(uprev).(extrapolation_scalars[n_curr + 1]) *
            sum(broadcast(*, T[1:(n_curr + 1)],
            eltype(uprev).(extrapolation_weights[1:(n_curr + 1),
                (n_curr + 1)]))) # Approximation of extrapolation order n_curr
    end

    # Save the latest approximation and update FSAL
    integrator.u = u
    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ImplicitEulerBarycentricExtrapolationCache)
    # cf. initialize! of MidpointCache
    @unpack k, fsalfirst = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

function perform_step!(integrator, cache::ImplicitEulerBarycentricExtrapolationCache,
        repeat_step = false)
    # Unpack all information needed
    @unpack t, uprev, dt, f, p = integrator
    alg = unwrap_alg(integrator, true)
    @unpack n_curr, u_temp1, u_temp2, utilde, res, T, fsalfirst, k, diff1, diff2 = cache
    @unpack u_temp3, u_temp4, k_tmps = cache
    # Coefficients for obtaining u
    @unpack extrapolation_weights, extrapolation_scalars = cache.coefficients
    # Coefficients for obtaining utilde
    @unpack extrapolation_weights_2, extrapolation_scalars_2 = cache.coefficients
    # Additional constant information
    @unpack subdividing_sequence = cache.coefficients
    @unpack sequence_factor = alg

    @unpack J, W, uf, tf, linsolve_tmps, jac_config = cache

    fill!(cache.Q, zero(eltype(cache.Q)))

    if integrator.opts.adaptive
        # Set up the order window
        # alg.min_order + 1 ≦ n_curr ≦ alg.max_order - 1 is enforced by step_*_controller!
        if !(alg.min_order + 1 <= n_curr <= alg.max_order - 1)
            error("Something went wrong while setting up the order window: $n_curr ∉ [$(alg.min_order+1),$(alg.max_order-1)].
            Please report this error  ")
        end
        win_min = n_curr - 1
        win_max = n_curr + 1

        # Set up the current extrapolation order
        cache.n_old = n_curr # Save the suggested order for step_*_controller!
        n_curr = win_min # Start with smallest order in the order window
    end

    #Compute the internal discretisations
    calc_J!(J, integrator, cache) # Store the calculated jac as it won't change in internal discretisation
    if !isthreaded(alg.threading)
        for i in 0:n_curr
            j_int = sequence_factor * subdividing_sequence[i + 1]
            dt_int = dt / j_int # Stepsize of the ith internal discretisation
            jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
            integrator.stats.nw += 1
            @.. broadcast=false u_temp2=uprev
            @.. broadcast=false linsolve_tmps[1]=fsalfirst

            linsolve = cache.linsolve[1]

            if !repeat_step
                linres = dolinsolve(integrator, linsolve; A = W[1],
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            else
                linres = dolinsolve(integrator, linsolve; A = nothing,
                    b = _vec(linsolve_tmps[1]), linu = _vec(k))
            end
            cache.linsolve[1] = linres.cache

            integrator.stats.nsolve += 1
            @.. broadcast=false u_temp1=u_temp2 - k # Euler starting step
            @.. broadcast=false diff1[1]=u_temp1 - u_temp2
            for j in 2:(j_int + 1)
                f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast=false linsolve_tmps[1]=k

                linsolve = cache.linsolve[1]
                if !repeat_step && j == 1
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                end
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false T[i + 1]=u_temp1 - k
                if (j == j_int + 1)
                    @.. broadcast=false T[i + 1]=0.25(T[i + 1] + 2 * u_temp1 + u_temp2)
                end
                @.. broadcast=false u_temp2=u_temp1
                @.. broadcast=false u_temp1=T[i + 1]
                if (i <= 1)
                    # Deuflhard Stability check for initial two sequences
                    @.. broadcast=false diff2[1]=u_temp1 - u_temp2
                    @.. broadcast=false diff2[1]=0.5 * (diff2[1] - diff1[1])
                    if (integrator.opts.internalnorm(diff1[1], t) <
                        integrator.opts.internalnorm(diff2[1], t))
                        # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                        integrator.force_stepfail = true
                        return
                    end
                end
                @.. broadcast=false diff1[1]=u_temp1 - u_temp2
            end
        end
    else
        if alg.sequence == :romberg
            # Compute solution by using maximum two threads for romberg sequence
            # One thread will fill T matrix till second last element and another thread will
            # fill last element of T matrix.
            # Romberg sequence --> 1, 2, 4, 8, ..., 2^(i)
            # 1 + 2 + 4 + ... + 2^(i-1) = 2^(i) - 1
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 1:2
                    startIndex = (i == 1) ? 0 : n_curr
                    endIndex = (i == 1) ? n_curr - 1 : n_curr

                    for index in startIndex:endIndex
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix,
                            dt_int_temp, J)
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=fsalfirst

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end
                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false k_tmps[Threads.threadid()]=-k_tmps[Threads.threadid()]
                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] -
                                                                        k_tmps[Threads.threadid()] # Euler starting step
                        @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                      u_temp4[Threads.threadid()]
                        for j in 2:(j_int_temp + 1)
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()]

                            linsolve = cache.linsolve[Threads.threadid()]

                            if !repeat_step && j == 1
                                linres = dolinsolve(integrator, linsolve;
                                    A = W[Threads.threadid()],
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            end
                            cache.linsolve[Threads.threadid()] = linres.cache

                            @.. broadcast=false T[index + 1]=u_temp3[Threads.threadid()] -
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            if (j == j_int_temp + 1)
                                @.. broadcast=false T[index + 1]=0.25(T[index + 1] +
                                                                      2 *
                                                                      u_temp3[Threads.threadid()] +
                                                                      u_temp4[Threads.threadid()])
                            end
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                              u_temp4[Threads.threadid()]
                                @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                              (diff2[Threads.threadid()] -
                                                                               diff1[Threads.threadid()])
                                if (integrator.opts.internalnorm(diff1[Threads.threadid()],
                                    t) <
                                    integrator.opts.internalnorm(diff2[Threads.threadid()],
                                    t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                          u_temp4[Threads.threadid()]
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        else
            let n_curr = n_curr, subdividing_sequence = subdividing_sequence, uprev = uprev,
                dt = dt, u_temp3 = u_temp3,
                u_temp4 = u_temp4, k_tmps = k_tmps, p = p, t = t, T = T

                @threaded alg.threading for i in 0:(n_curr ÷ 2)
                    indices = i != n_curr - i ? (i, n_curr - i) : (-1, n_curr - i)
                    for index in indices
                        index == -1 && continue
                        j_int_temp = sequence_factor * subdividing_sequence[index + 1]
                        dt_int_temp = dt / j_int_temp # Stepsize of the ith internal discretisation
                        jacobian2W!(W[Threads.threadid()], integrator.f.mass_matrix,
                            dt_int_temp, J)
                        @.. broadcast=false u_temp4[Threads.threadid()]=uprev
                        @.. broadcast=false linsolve_tmps[Threads.threadid()]=fsalfirst

                        linsolve = cache.linsolve[Threads.threadid()]

                        if !repeat_step
                            linres = dolinsolve(integrator, linsolve;
                                A = W[Threads.threadid()],
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        else
                            linres = dolinsolve(integrator, linsolve; A = nothing,
                                b = _vec(linsolve_tmps[Threads.threadid()]),
                                linu = _vec(k_tmps[Threads.threadid()]))
                        end
                        cache.linsolve[Threads.threadid()] = linres.cache

                        @.. broadcast=false u_temp3[Threads.threadid()]=u_temp4[Threads.threadid()] -
                                                                        k_tmps[Threads.threadid()] # Euler starting step
                        @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                      u_temp4[Threads.threadid()]
                        for j in 2:(j_int_temp + 1)
                            f(k_tmps[Threads.threadid()],
                                cache.u_temp3[Threads.threadid()],
                                p, t + (j - 1) * dt_int_temp)
                            @.. broadcast=false linsolve_tmps[Threads.threadid()]=k_tmps[Threads.threadid()]

                            linsolve = cache.linsolve[Threads.threadid()]

                            if (!repeat_step && j == 1)
                                linres = dolinsolve(integrator, linsolve;
                                    A = W[Threads.threadid()],
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            else
                                linres = dolinsolve(integrator, linsolve; A = nothing,
                                    b = _vec(linsolve_tmps[Threads.threadid()]),
                                    linu = _vec(k_tmps[Threads.threadid()]))
                            end
                            cache.linsolve[Threads.threadid()] = linres.cache

                            @.. broadcast=false T[index + 1]=u_temp3[Threads.threadid()] -
                                                             k_tmps[Threads.threadid()] # Explicit Midpoint rule
                            if (j == j_int_temp + 1)
                                @.. broadcast=false T[index + 1]=0.25(T[index + 1] +
                                                                      2 *
                                                                      u_temp3[Threads.threadid()] +
                                                                      u_temp4[Threads.threadid()])
                            end
                            @.. broadcast=false u_temp4[Threads.threadid()]=u_temp3[Threads.threadid()]
                            @.. broadcast=false u_temp3[Threads.threadid()]=T[index + 1]
                            if (index <= 1)
                                # Deuflhard Stability check for initial two sequences
                                @.. broadcast=false diff2[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                              u_temp4[Threads.threadid()]
                                @.. broadcast=false diff2[Threads.threadid()]=0.5 *
                                                                              (diff2[Threads.threadid()] -
                                                                               diff1[Threads.threadid()])
                                if (integrator.opts.internalnorm(diff1[Threads.threadid()],
                                    t) <
                                    integrator.opts.internalnorm(diff2[Threads.threadid()],
                                    t))
                                    # Divergence of iteration, overflow is possible. Force fail and start with smaller step
                                    integrator.force_stepfail = true
                                    return
                                end
                            end
                            @.. broadcast=false diff1[Threads.threadid()]=u_temp3[Threads.threadid()] -
                                                                          u_temp4[Threads.threadid()]
                        end
                    end
                    integrator.force_stepfail ? break : continue
                end
            end
        end
    end

    if integrator.force_stepfail
        return
    end

    if integrator.opts.adaptive
        # Compute all information relating to an extrapolation order ≦ win_min
        for i in (win_min - 1):win_min

            #integrator.u .= extrapolation_scalars[i+1] * sum( broadcast(*, cache.T[1:(i+1)], extrapolation_weights[1:(i+1), (i+1)]) ) # Approximation of extrapolation order i
            #cache.utilde .= extrapolation_scalars_2[i] * sum( broadcast(*, cache.T[2:(i+1)], extrapolation_weights_2[1:i, i]) ) # and its internal counterpart

            u_temp1 .= false
            u_temp2 .= false
            for j in 1:(i + 1)
                @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (i + 1)]
            end
            for j in 2:(i + 1)
                @.. broadcast=false u_temp2+=cache.T[j] * extrapolation_weights_2[j - 1, i]
            end
            @.. broadcast=false integrator.u=extrapolation_scalars[i + 1] * u_temp1
            @.. broadcast=false cache.utilde=extrapolation_scalars_2[i] * u_temp2

            calculate_residuals!(cache.res, integrator.u, cache.utilde,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(cache.res, t)
            cache.n_curr = i # Update cache's n_curr for stepsize_controller_internal!
            stepsize_controller_internal!(integrator, alg) # Update cache.Q
        end

        # Check if an approximation of some order in the order window can be accepted
        # Make sure a stepsize scaling factor of order (alg.min_order + 1) is provided for the step_*_controller!
        while n_curr <= win_max
            EEst1 = one(integrator.EEst)
            for i in (n_curr + 2):(win_max + 1)
                EEst1 *= subdividing_sequence[i] / subdividing_sequence[1]
            end
            EEst1 *= EEst1

            #@show integrator.opts.internalnorm(integrator.u - cache.utilde,t)
            if accept_step_controller(integrator, integrator.opts.controller)
                # Accept current approximation u of order n_curr
                break
            elseif (n_curr < alg.min_order + 1) || integrator.EEst <= EEst1
                # Reject current approximation order but pass convergence monitor
                # Compute approximation of order (n_curr + 1)
                n_curr = n_curr + 1
                cache.n_curr = n_curr

                # Update cache.T
                j_int = sequence_factor * subdividing_sequence[n_curr + 1]
                dt_int = dt / j_int # Stepsize of the new internal discretisation
                jacobian2W!(W[1], integrator.f.mass_matrix, dt_int, J)
                integrator.stats.nw += 1
                @.. broadcast=false u_temp2=uprev
                @.. broadcast=false linsolve_tmps[1]=fsalfirst

                linsolve = cache.linsolve[1]

                if !repeat_step
                    linres = dolinsolve(integrator, linsolve; A = W[1],
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                else
                    linres = dolinsolve(integrator, linsolve; A = nothing,
                        b = _vec(linsolve_tmps[1]), linu = _vec(k))
                end
                cache.linsolve[1] = linres.cache

                integrator.stats.nsolve += 1
                @.. broadcast=false k=-k
                @.. broadcast=false u_temp1=u_temp2 + k # Euler starting step
                for j in 2:(j_int + 1)
                    f(k, cache.u_temp1, p, t + (j - 1) * dt_int)
                    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                    @.. broadcast=false linsolve_tmps[1]=k

                    linsolve = cache.linsolve[1]
                    linres = dolinsolve(integrator, linsolve; b = _vec(linsolve_tmps[1]),
                        linu = _vec(k))
                    cache.linsolve[1] = linres.cache

                    integrator.stats.nsolve += 1
                    @.. broadcast=false T[n_curr + 1]=u_temp1 - k # Explicit Midpoint rule
                    if (j == j_int + 1)
                        @.. broadcast=false T[n_curr + 1]=0.25(T[n_curr + 1] + 2 * u_temp1 +
                                                               u_temp2)
                    end
                    @.. broadcast=false u_temp2=u_temp1
                    @.. broadcast=false u_temp1=T[n_curr + 1]
                end

                # Update u, integrator.EEst and cache.Q
                #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
                #cache.utilde .= extrapolation_scalars_2[n_curr] * sum( broadcast(*, cache.T[2:(n_curr+1)], extrapolation_weights_2[1:n_curr, n_curr]) ) # and its internal counterpart

                u_temp1 .= false
                u_temp2 .= false
                for j in 1:(n_curr + 1)
                    @.. broadcast=false u_temp1+=cache.T[j] *
                                                 extrapolation_weights[j, (n_curr + 1)]
                end
                for j in 2:(n_curr + 1)
                    @.. broadcast=false u_temp2+=cache.T[j] *
                                                 extrapolation_weights_2[j - 1, n_curr]
                end
                @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
                @.. broadcast=false cache.utilde=extrapolation_scalars_2[n_curr] * u_temp2

                calculate_residuals!(cache.res, integrator.u, cache.utilde,
                    integrator.opts.abstol, integrator.opts.reltol,
                    integrator.opts.internalnorm, t)
                integrator.EEst = integrator.opts.internalnorm(cache.res, t)
                stepsize_controller_internal!(integrator, alg) # Update cache.Q
            else
                # Reject the current approximation and not pass convergence monitor
                break
            end
        end
    else

        #integrator.u .= extrapolation_scalars[n_curr+1] * sum( broadcast(*, cache.T[1:(n_curr+1)], extrapolation_weights[1:(n_curr+1), (n_curr+1)]) ) # Approximation of extrapolation order n_curr
        u_temp1 .= false
        for j in 1:(n_curr + 1)
            @.. broadcast=false u_temp1+=cache.T[j] * extrapolation_weights[j, (n_curr + 1)]
        end
        @.. broadcast=false integrator.u=extrapolation_scalars[n_curr + 1] * u_temp1
    end

    f(cache.k, integrator.u, p, t + dt) # Update FSAL
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

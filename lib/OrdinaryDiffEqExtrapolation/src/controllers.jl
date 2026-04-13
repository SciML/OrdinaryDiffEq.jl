@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        # Extrapolation methods
        mutable struct ExtrapolationController{QT} <: AbstractLegacyController
            beta1::QT
        end

        struct NewExtrapolationController{QT} <: AbstractController
            qmin::QT
            qmax::QT
        end
    end
else
    @eval begin
        mutable struct ExtrapolationController{QT} <: AbstractController
            beta1::QT
        end
    end
end

function reset_alg_dependent_opts!(controller::ExtrapolationController, alg1, alg2)
    if controller.beta1 == beta1_default(alg1, beta2_default(alg1))
        controller.beta1 = beta1_default(alg2, beta2_default(alg2))
    end
    return nothing
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function NewExtrapolationController(QT, alg; qmin = nothing, qmax = nothing, gamma = nothing)
            return NewExtrapolationController(
                QT(qmin === nothing ? qmin_default(alg) : qmin),
                QT(qmax === nothing ? qmax_default(alg) : qmax),
            )
        end

        mutable struct ExtrapolationControllerCache{QT, UT} <: AbstractControllerCache
            controller::NewExtrapolationController{QT}
            beta1::QT
            gamma::QT
            atmp::UT
        end

        function setup_controller_cache(alg, atmp, controller::NewExtrapolationController{T}) where {T}
            return ExtrapolationControllerCache(
                controller,
                T(1),
                T(1),
                atmp,
            )
        end
    end
end

# TODO replace this when done - right now this holds the controller cache!
stepsize_controller_internal!(integrator, alg) = stepsize_controller_internal!(integrator, integrator.opts.controller, alg)
stepsize_predictor!(integrator, alg, n_new) = stepsize_predictor!(integrator, integrator.opts.controller, alg, n_new)
# stepsize_controller_internal!(integrator, alg) = stepsize_controller_internal!(integrator, integrator.cache.controller_cache, alg)
# stepsize_predictor!(integrator, alg, n_new) = stepsize_predictor!(integrator, integrator.cache.controller_cache, alg, n_new)

@inline function stepsize_controller!(
        integrator,
        controller::ExtrapolationController,
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
        }
    )
    # Dummy function
    # ExtrapolationMidpointDeuflhard's stepsize scaling is stored in the cache;
    # it is computed by  stepsize_controller_internal! (in perform_step!) resp. stepsize_predictor!
    # (in step_accept_controller! and step_reject_controller!)
    return zero(typeof(integrator.opts.qmax))
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        @inline function stepsize_controller!(
                integrator,
                cache::ExtrapolationControllerCache,
                alg::Union{
                    ExtrapolationMidpointDeuflhard,
                    ImplicitDeuflhardExtrapolation,
                }
            )
            # Dummy function
            # ExtrapolationMidpointDeuflhard's stepsize scaling is stored in the cache;
            # it is computed by  stepsize_controller_internal! (in perform_step!) resp. stepsize_predictor!
            # (in step_accept_controller! and step_reject_controller!)
            return zero(typeof(cache.controller.qmax))
        end
    end
end

function stepsize_controller_internal!(
        integrator,
        controller::ExtrapolationController,
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
        }
    )
    # Standard step size controller
    # Compute and save the stepsize scaling based on the latest error estimate of the current order
    (; controller) = integrator.opts

    if iszero(integrator.EEst)
        q = inv(integrator.opts.qmax)
    else
        # Update gamma and beta1
        controller.beta1 = typeof(controller.beta1)(1 // (2integrator.cache.n_curr + 1))
        integrator.opts.gamma = FastPower.fastpower(
            typeof(integrator.opts.gamma)(1 // 4),
            controller.beta1
        )
        # Compute new stepsize scaling
        qtmp = FastPower.fastpower(integrator.EEst, controller.beta1) /
            integrator.opts.gamma
        @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
    end
    return integrator.cache.Q[integrator.cache.n_curr - alg.min_order + 1] = q
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function stepsize_controller_internal!(
                integrator,
                cache::ExtrapolationControllerCache,
                alg::Union{
                    ExtrapolationMidpointDeuflhard,
                    ImplicitDeuflhardExtrapolation,
                }
            )
            # Standard step size controller
            # Compute and save the stepsize scaling based on the latest error estimate of the current order
            (; controller) = cache

            if iszero(integrator.EEst)
                q = inv(controller.qmax)
            else
                # Update gamma and beta1
                cache.beta1 = typeof(cache.beta1)(1 // (2integrator.cache.n_curr + 1))
                cache.gamma = FastPower.fastpower(
                    typeof(cache.gamma)(1 // 4),
                    cache.beta1
                )
                # Compute new stepsize scaling
                qtmp = FastPower.fastpower(integrator.EEst, cache.beta1) /
                    cache.gamma
                @fastmath q = max(inv(controller.qmax), min(inv(controller.qmin), qtmp))
            end
            return integrator.cache.Q[integrator.cache.n_curr - alg.min_order + 1] = q
        end
    end
end

function stepsize_predictor!(
        integrator,
        controller::ExtrapolationController,
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
        }, n_new::Int
    )
    # Compute and save the stepsize scaling for order n_new based on the latest error estimate of the current order.
    (; controller) = integrator.opts

    if iszero(integrator.EEst)
        q = inv(integrator.opts.qmax)
    else
        # Initialize
        (; t, EEst) = integrator
        (; stage_number) = integrator.cache
        tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Deuflhard's approach relies on EEstD ≈ ||relTol||
        s_curr = stage_number[integrator.cache.n_curr - alg.min_order + 1]
        s_new = stage_number[n_new - alg.min_order + 1]
        # Update gamma and beta1
        controller.beta1 = typeof(controller.beta1)(1 // (2integrator.cache.n_curr + 1))
        integrator.opts.gamma = FastPower.fastpower(
            typeof(integrator.opts.gamma)(1 // 4),
            controller.beta1
        )
        # Compute new stepsize scaling
        qtmp = EEst *
            FastPower.fastpower(
            FastPower.fastpower(tol, (1.0 - s_curr / s_new)),
            controller.beta1
        ) / integrator.opts.gamma
        @fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
    end
    return integrator.cache.Q[n_new - alg.min_order + 1] = q
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function stepsize_predictor!(
                integrator,
                cache::ExtrapolationControllerCache,
                alg::Union{
                    ExtrapolationMidpointDeuflhard,
                    ImplicitDeuflhardExtrapolation,
                }, n_new::Int
            )
            # Compute and save the stepsize scaling for order n_new based on the latest error estimate of the current order.
            (; controller) = cache

            if iszero(integrator.EEst)
                q = inv(controller.qmax)
            else
                # Initialize
                (; t, EEst) = integrator
                (; stage_number) = integrator.cache
                tol = integrator.opts.internalnorm(integrator.opts.reltol, t) # Deuflhard's approach relies on EEstD ≈ ||relTol||
                s_curr = stage_number[integrator.cache.n_curr - alg.min_order + 1]
                s_new = stage_number[n_new - alg.min_order + 1]
                # Update gamma and beta1
                cache.beta1 = typeof(cache.beta1)(1 // (2integrator.cache.n_curr + 1))
                cache.gamma = FastPower.fastpower(
                    typeof(cache.gamma)(1 // 4),
                    cache.beta1
                )
                # Compute new stepsize scaling
                qtmp = EEst *
                    FastPower.fastpower(
                    FastPower.fastpower(tol, (1.0 - s_curr / s_new)),
                    cache.beta1
                ) / cache.gamma
                @fastmath q = max(inv(controller.qmax), min(inv(controller.qmin), qtmp))
            end
            return integrator.cache.Q[n_new - alg.min_order + 1] = q
        end
    end
end

function step_accept_controller!(
        integrator,
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
        }, q
    )
    # Compute new order and stepsize, return new stepsize
    (; min_order, max_order) = alg
    (; n_curr, n_old, Q) = integrator.cache
    s = integrator.cache.stage_number

    # Compute new order based on available quantities
    tmp = (min_order:n_curr) .- min_order .+ 1 # Index range of quantities computed so far
    dt_new = Vector{eltype(Q)}(undef, length(tmp) + 1)
    dt_new[1:(end - 1)] = integrator.dt ./ Q[tmp] # Store for the possible new stepsizes
    dtmin = timedepentdtmin(integrator)
    dt_new[1:(end - 1)] = max.(
        dtmin,
        min.(abs(integrator.opts.dtmax), abs.(dt_new[1:(end - 1)]))
    ) # Safety scaling

    # n_new is the most efficient order of the last step
    work = s[tmp] ./ dt_new[1:(end - 1)]
    n_new = argmin(work) + min_order - 1

    # Check if n_new may be increased
    if n_new == n_curr < min(max_order, n_old + 1) # cf. win_max in perform_step! of the last step
        # Predict stepsize scaling for order (n_new + 1)
        stepsize_predictor!(integrator, alg, n_new + 1) # Update cache.Q

        # Compute and scale the corresponding stepsize
        dt_new[end] = integrator.dt ./ Q[tmp[end] + 1]
        dt_new[end] = max(dtmin, min(abs(integrator.opts.dtmax), abs.(dt_new[end])))

        # Check if (n_new  + 1) would have been more efficient than n_new
        if work[end] > s[tmp[end] + 1] / dt_new[end]
            n_new = n_new + 1
        end
    end

    integrator.cache.n_curr = n_new
    return dt_new[n_new - min_order + 1]
end

function step_reject_controller!(
        integrator,
        alg::Union{
            ExtrapolationMidpointDeuflhard,
            ImplicitDeuflhardExtrapolation,
        }
    )
    # Compute and save reduced stepsize dt_red of order n_old
    # Use the latest error estimate to predict dt_red if an estimate of order n_old is not available
    if integrator.cache.n_curr < integrator.cache.n_old
        stepsize_predictor!(integrator, alg, integrator.cache.n_old) # Update cache.Q
    end
    integrator.cache.n_curr = integrator.cache.n_old # Reset order for redoing the rejected step
    dt_red = integrator.dt /
        integrator.cache.Q[integrator.cache.n_old - integrator.alg.min_order + 1]
    dtmin = timedepentdtmin(integrator)
    dt_red = integrator.tdir * max(dtmin, min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
    return integrator.dt = dt_red
end

@inline function stepsize_controller!(
        integrator,
        alg::Union{
            ExtrapolationMidpointHairerWanner,
            ImplicitHairerWannerExtrapolation,
            ImplicitEulerExtrapolation,
            ImplicitEulerBarycentricExtrapolation,
        }
    )
    # Dummy function
    # ExtrapolationMidpointHairerWanner's stepsize scaling is stored in the cache;
    # it is computed by  stepsize_controller_internal! (in perform_step!), step_accept_controller! or step_reject_controller!
    return zero(typeof(integrator.opts.qmax))
end

function stepsize_controller_internal!(
        integrator,
        controller::ExtrapolationController,
        alg::Union{
            ExtrapolationMidpointHairerWanner,
            ImplicitHairerWannerExtrapolation,
            ImplicitEulerExtrapolation,
            ImplicitEulerBarycentricExtrapolation,
        }
    )
    # Standard step size controller
    # Compute and save the stepsize scaling based on the latest error estimate of the current order
    (; controller) = integrator.opts

    return if alg isa
            Union{
            ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation,
            ImplicitHairerWannerExtrapolation,
        }
        if iszero(integrator.EEst)
            q = inv(integrator.opts.qmax)
        else
            # Update gamma and beta1
            if alg isa ImplicitHairerWannerExtrapolation
                controller.beta1 = typeof(controller.beta1)(
                    1 //
                        (2integrator.cache.n_curr + 1)
                )
            elseif alg isa ImplicitEulerExtrapolation
                controller.beta1 = typeof(controller.beta1)(1 // (integrator.cache.n_curr))
            else
                controller.beta1 = typeof(controller.beta1)(
                    1 //
                        (integrator.cache.n_curr - 1)
                )
            end
            integrator.opts.gamma = FastPower.fastpower(
                typeof(integrator.opts.gamma)(
                    65 //
                        100
                ),
                controller.beta1
            )
            # Compute new stepsize scaling
            qtmp = FastPower.fastpower(integrator.EEst, controller.beta1) /
                (integrator.opts.gamma)
            @fastmath q = max(
                inv(integrator.opts.qmax),
                min(inv(integrator.opts.qmin), qtmp)
            )
        end
        integrator.cache.Q[integrator.cache.n_curr + 1] = q
    else
        if iszero(integrator.EEst)
            q = inv(integrator.opts.qmax)
        else
            # Update gamma and beta1
            controller.beta1 = typeof(controller.beta1)(1 // (2integrator.cache.n_curr + 1))
            integrator.opts.gamma = FastPower.fastpower(
                typeof(integrator.opts.gamma)(
                    65 //
                        100
                ),
                controller.beta1
            )
            # Compute new stepsize scaling
            qtmp = FastPower.fastpower(integrator.EEst, controller.beta1) /
                integrator.opts.gamma
            @fastmath q = max(
                inv(integrator.opts.qmax),
                min(inv(integrator.opts.qmin), qtmp)
            )
        end
        integrator.cache.Q[integrator.cache.n_curr + 1] = q
    end
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function stepsize_controller_internal!(
                integrator,
                cache::ExtrapolationControllerCache,
                alg::Union{
                    ExtrapolationMidpointHairerWanner,
                    ImplicitHairerWannerExtrapolation,
                    ImplicitEulerExtrapolation,
                    ImplicitEulerBarycentricExtrapolation,
                }
            )
            # Standard step size controller
            # Compute and save the stepsize scaling based on the latest error estimate of the current order
            (; controller) = cache

            return if alg isa
                    Union{
                    ImplicitEulerExtrapolation, ImplicitEulerBarycentricExtrapolation,
                    ImplicitHairerWannerExtrapolation,
                }
                if iszero(integrator.EEst)
                    q = inv(controller.qmax)
                else
                    # Update gamma and beta1
                    if alg isa ImplicitHairerWannerExtrapolation
                        cache.beta1 = typeof(cache.beta1)(
                            1 //
                                (2integrator.cache.n_curr + 1)
                        )
                    elseif alg isa ImplicitEulerExtrapolation
                        cache.beta1 = typeof(cache.beta1)(1 // (integrator.cache.n_curr))
                    else
                        cache.beta1 = typeof(cache.beta1)(
                            1 //
                                (integrator.cache.n_curr - 1)
                        )
                    end
                    cache.gamma = FastPower.fastpower(
                        typeof(cache.gamma)(
                            65 //
                                100
                        ),
                        cache.beta1
                    )
                    # Compute new stepsize scaling
                    qtmp = FastPower.fastpower(integrator.EEst, cache.beta1) /
                        (cache.gamma)
                    @fastmath q = max(
                        inv(controller.qmax),
                        min(inv(controller.qmin), qtmp)
                    )
                end
                integrator.cache.Q[integrator.cache.n_curr + 1] = q
            else
                if iszero(integrator.EEst)
                    q = inv(controller.qmax)
                else
                    # Update gamma and beta1
                    cache.beta1 = typeof(cache.beta1)(1 // (2integrator.cache.n_curr + 1))
                    cache.gamma = FastPower.fastpower(
                        typeof(cache.gamma)(
                            65 //
                                100
                        ),
                        cache.beta1
                    )
                    # Compute new stepsize scaling
                    qtmp = FastPower.fastpower(integrator.EEst, cache.beta1) /
                        cache.gamma
                    @fastmath q = max(
                        inv(controller.qmax),
                        min(inv(controller.qmin), qtmp)
                    )
                end
                integrator.cache.Q[integrator.cache.n_curr + 1] = q
            end
        end
    end
end

function step_accept_controller!(
        integrator,
        alg::Union{
            ExtrapolationMidpointHairerWanner,
            ImplicitHairerWannerExtrapolation,
            ImplicitEulerExtrapolation,
            ImplicitEulerBarycentricExtrapolation,
        }, q
    )
    # Compute new order and stepsize, return new stepsize
    (; min_order, max_order) = alg
    (; n_curr, n_old, Q, sigma, work, dt_new) = integrator.cache
    s = integrator.cache.stage_number

    # Compute new order based on available quantities
    win_min_old = min(n_old, n_curr) - 1 # cf. win_min in perform_step! of the last step
    tmp = win_min_old:(max(n_curr, n_old) + 1) # Index range for the new order
    fill!(dt_new, zero(eltype(dt_new)))
    @.. broadcast = false Q = integrator.dt / Q
    copyto!(dt_new, win_min_old, Q, win_min_old, (max(n_curr, n_old) + 1) - win_min_old + 1)
    @.. broadcast = false Q = integrator.dt / Q
    dtmin = timedepentdtmin(integrator)
    fill!(work, zero(eltype(work))) # work[n] is the work for order (n-1)
    for i in tmp
        work[i] = s[i] / dt_new[i]
    end
    # Order selection
    n_new = n_old
    if n_curr == min_order # Enforce min_order + 1 ≦ n_new
        n_new = min_order + 1
    else
        if n_curr <= n_old
            if work[n_curr - 1] < sigma * work[n_curr]
                n_new = max(n_curr - 1, n_old - 1, min_order + 1) # Enforce min_order + 1≦ n_new
            elseif work[n_curr] < sigma * work[n_curr - 1]
                n_new = min(n_curr + 1, max_order - 1) # Enforce n_new ≦ max_order - 1
            else
                n_new = n_curr # min_order + 1 ≦ n_curr
            end
        else
            if work[n_old] < sigma * work[n_old + 1]
                n_new = max(n_old - 1, min_order + 1)  # Enforce min_order + 1 ≦ n_new
            end
            if work[n_curr + 1] < sigma * work[n_new + 1]
                n_new = min(n_new + 1, max_order - 1) # Enforce n_new ≦ max_order - 1
            end
        end
    end
    integrator.cache.n_curr = n_new

    # Stepsize selection
    if n_new == n_curr + 1
        # Compute the new stepsize of order n_new based on the optimal stepsize of order n_curr
        dt_new[n_new + 1] = s[n_curr + 2] / s[n_curr + 1] * dt_new[n_curr + 1]
        dt_new[n_new + 1] = max(
            dtmin,
            min(abs(integrator.opts.dtmax), abs(dt_new[n_new + 1]))
        )
    end
    return dt_new[n_new + 1]
end

function step_reject_controller!(
        integrator,
        alg::Union{
            ExtrapolationMidpointHairerWanner,
            ImplicitHairerWannerExtrapolation,
            ImplicitEulerExtrapolation,
            ImplicitEulerBarycentricExtrapolation,
        }
    )
    # Compute and save order and stepsize for redoing the current step
    (; n_old, n_curr, Q) = integrator.cache

    # Order selection
    n_red = n_old
    if n_curr == n_old - 1
        n_red = max(alg.min_order + 1, n_old - 1) # Enforce min_order + 1 ≦ n_red
    end
    integrator.cache.n_curr = n_red

    # Stepsize selection
    dt_red = integrator.dt / Q[n_red + 1]
    dtmin = timedepentdtmin(integrator)
    dt_red = integrator.tdir * max(dtmin, min(abs(integrator.opts.dtmax), abs(dt_red))) # Safety scaling
    return integrator.dt = dt_red
end

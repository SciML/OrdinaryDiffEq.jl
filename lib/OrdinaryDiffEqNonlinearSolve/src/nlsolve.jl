@inline eps_around_one(θ::T) where {T} = 100sqrt(eps(one(θ)))

# A globalized inner solver (TrustRegion and friends) can reject its trial step: `step!`
# returns with the iterate exactly unmoved and the cache not terminated, having only shrunk
# its trust region for the next attempt. The zero displacement such a `step!` leaves behind
# is not convergence evidence — the convergence tests below would read it as a perfect solve
# (`ndz < 1e-5` on the first iteration, `η·ndz = 0 < κ` on later ones) and accept the stage
# with no correction applied (#3817). A cache that instead *terminated* at zero displacement
# either converged exactly (which genuinely is convergence) or failed, and failures were
# already turned into `Inf` by `compute_step!`, so termination is the discriminator between
# "already at the root" and "has not decided anything yet".
_rejected_trial_step(nlsolver, ndz) = false
function _rejected_trial_step(nlsolver::NLSolver{<:NonlinearSolveAlg}, ndz)
    nlcache = nlsolver.cache.cache
    # A no-init cache is `solve!`-driven, so it has no trial-step state to reject (and no
    # `force_stop`/`nsteps` for `not_terminated` to read): a zero `ndz` there means a
    # complete inner solve returned the iterate unchanged, which is genuine convergence.
    nlcache isa NonlinearSolveNoInitCache && return false
    return iszero(ndz) && NonlinearSolveBase.not_terminated(nlcache)
end

"""
    nlsolve!(nlsolver::AbstractNLSolver, integrator)

Solve

```math
dt⋅f(innertmp + γ⋅z, p, t + c⋅dt) + outertmp = z
```

where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.

Whether `innertmp` and `outertmp` is used for the evaluation is controlled by setting `nlsolver.method`.
In both cases the variable name is actually `nlsolver.tmp`.
"""
function nlsolve!(
        nlsolver::NL, integrator::SciMLBase.DEIntegrator,
        cache = nothing, repeat_step = false
    ) where {NL <: AbstractNLSolver}
    always_new = is_always_new(nlsolver)
    check_div′ = check_div(nlsolver)
    @label REDO
    # Initialize γW for JET
    γW = one(integrator.dt)
    if isnewton(nlsolver)
        cache === nothing &&
            throw(ArgumentError("cache is not passed to `nlsolve!` when using NLNewton"))
        if nlsolver.method === DIRK
            γW = nlsolver.γ * integrator.dt
        else
            γW = nlsolver.γ * integrator.dt / nlsolver.α
        end
        always_new || update_W!(nlsolver, integrator, cache, γW, repeat_step)
    end

    (; maxiters, κ, fast_convergence_cutoff) = nlsolver

    initialize!(nlsolver, integrator)
    nlsolver.status = check_div′ ? Divergence : Convergence
    η = get_new_W!(nlsolver) ? initial_η(nlsolver, integrator) : nlsolver.ηold

    ndz = one(η)
    for iter in 1:maxiters
        if always_new && isnewton(nlsolver)
            if ArrayInterface.ismutable(integrator.u)
                @.. integrator.u = integrator.uprev + nlsolver.γ * nlsolver.z
            else
                integrator.u = @.. integrator.uprev + nlsolver.γ * nlsolver.z
            end
            update_W!(nlsolver, integrator, cache, γW, repeat_step, (true, true))
        end
        nlsolver.iter = iter

        # compute next step and calculate norm of residuals
        ndzprev = ndz
        if isnewton(nlsolver)
            # Newton solve requires γW in order to update W
            ndz = compute_step!(nlsolver, integrator, γW)
        else
            ndz = compute_step!(nlsolver, integrator)
        end
        if !isfinite(ndz)
            @SciMLMessage(
                lazy"Newton iteration diverged: residual norm is not finite (ndz = $(ndz))",
                integrator.opts.verbose, :newton_convergence
            )
            nlsolver.status = Divergence
            nlsolver.nfails += 1
            break
        end

        if _rejected_trial_step(nlsolver, ndz)
            @SciMLMessage(
                lazy"Inner nonlinear solver rejected its trial step (iter = $(iter)); the unmoved iterate is not treated as convergence",
                integrator.opts.verbose, :newton_convergence
            )
            # No new iterate exists to judge: skip the convergence and divergence
            # bookkeeping (a zero `ndz` would also poison the next iteration's `θ`)
            # and let the inner solver retry with its shrunken trust region. If it
            # never moves, the loop runs out and the step fails as unconverged.
            ndz = ndzprev
            continue
        end

        has_prev_θ = hasfield(NL, :prev_θ)
        prev_θ = has_prev_θ ? nlsolver.prev_θ : one(ndz)

        # check divergence (not in initial step)
        if iter > 1
            θ = prev_θ = has_prev_θ ? max(0.3 * prev_θ, ndz / ndzprev) : ndz / ndzprev

            # When one Newton iteration basically does nothing, it's likely that we
            # are at the precision limit of floating point number. Thus, we just call
            # it convergence/divergence according to `ndz` directly.
            if abs(θ - one(θ)) <= eps_around_one(θ)
                if ndz <= one(ndz)
                    @SciMLMessage(
                        lazy"Newton iteration converged at floating point limit: θ ≈ 1.0, ndz = $(ndz)",
                        integrator.opts.verbose, :convergence_limit
                    )
                    nlsolver.status = Convergence
                    nlsolver.nfails = 0
                    break
                elseif check_div′
                    nlsolver.status = Divergence
                    nlsolver.nfails += 1
                    break
                end
            end

            # divergence
            if check_div′ && θ > 2
                @SciMLMessage(
                    lazy"Newton iteration diverging: θ = $(θ) > 2, ndz = $(ndz), ndzprev = $(ndzprev)",
                    integrator.opts.verbose, :newton_convergence
                )
                nlsolver.status = Divergence
                nlsolver.nfails += 1
                break
            end
        else
            if !integrator.accept_step
                prev_θ = one(prev_θ)
            end
            θ = prev_θ
        end

        if has_prev_θ
            nlsolver.prev_θ = prev_θ
        end

        apply_step!(nlsolver, integrator)

        # check for convergence
        η = SciMLBase.value(θ / (1 - θ))
        # don't trust θ for non-adaptive on first iter because the solver doesn't provide feedback
        # for us to know whether our previous nlsolve converged sufficiently well
        check_η_convergence = (
            iter > 1 ||
                (
                isadaptive(integrator.alg) &&
                    (
                    isnewton(nlsolver) ||
                        (
                        nlsolver.alg isa NonlinearSolveAlg &&
                            nlsolver.cache.W !== nothing
                    )
                )
            )
        )
        if (iter == 1 && ndz < 1.0e-5) ||
                (check_η_convergence && η >= zero(η) && η * ndz < κ)
            @SciMLMessage(
                lazy"Newton iteration converged in $(iter) iterations: η = $(η), ndz = $(ndz)",
                integrator.opts.verbose, :newton_iterations
            )
            nlsolver.status = Convergence
            nlsolver.nfails = 0
            break
        end
    end

    if isnewton(nlsolver) && nlsolver.status == Divergence &&
            !isJcurrent(nlsolver, integrator)
        @SciMLMessage(
            lazy"Newton iteration failed with stale Jacobian, retrying with fresh Jacobian",
            integrator.opts.verbose, :newton_convergence
        )
        nlsolver.status = TryAgain
        nlsolver.nfails += 1
        always_new || @goto REDO
    end

    nlsolver.ηold = η
    return postamble!(nlsolver, integrator)
end

## default implementations

initialize!(::AbstractNLSolver, integrator::SciMLBase.DEIntegrator) = nothing

"""
    initial_η(nlsolver, integrator) -> η

Return the initial convergence-rate estimate `η` for a fresh nonlinear solve.
Functional/Anderson solvers reuse the previous `ηold`; the Newton solver method
derives it from the tolerance. Consumed by [`nlsolve!`](@ref).
"""
function initial_η(nlsolver::NLSolver, integrator)
    return max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
end

function apply_step!(
        nlsolver::NLSolver{algType, iip},
        integrator::SciMLBase.DEIntegrator
    ) where {algType, iip}
    if iip
        @.. broadcast = false nlsolver.z = nlsolver.ztmp
    else
        nlsolver.z = nlsolver.ztmp
    end

    return nothing
end

function SciMLBase.postamble!(nlsolver::NLSolver, integrator::SciMLBase.DEIntegrator)
    if SciMLBase.has_stats(integrator)
        integrator.stats.nnonliniter += nlsolver.iter

        if nlsolvefail(nlsolver)
            integrator.stats.nnonlinconvfail += 1
        end
    end
    integrator.force_stepfail = nlsolvefail(nlsolver)
    setfirststage!(nlsolver, false)
    isnewton(nlsolver) && (nlsolver.cache.firstcall = false)

    return nlsolver.z
end

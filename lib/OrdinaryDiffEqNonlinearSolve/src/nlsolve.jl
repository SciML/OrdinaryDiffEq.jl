@inline eps_around_one(θ::T) where {T} = 100sqrt(eps(one(θ)))

# A globalized inner solver can return from `step!` with the iterate all but unmoved: a
# TrustRegion that rejected its trial step and only shrank its radius (#3817), a line search
# that collapsed to a near-zero step length. The convergence tests below would read that
# displacement as a perfect solve (`ndz < 1e-5` on the first iteration, `η·ndz ≈ 0 < κ` on
# later ones) and accept the stage with no correction applied.
#
# Two disjoint symptoms. An exactly unmoved iterate from a cache that has not terminated has
# decided nothing yet (a cache that terminated at zero displacement either reached the root
# exactly or failed, and `compute_step!` turns failures into `Inf`). A displacement that is
# merely below roundoff is not self-evidently either, so `compute_step!` puts the question to
# the residual — see `stalled_inner_step`.
_uninformative_step(nlsolver, ndz) = false
function _uninformative_step(nlsolver::NLSolver{<:NonlinearSolveAlg}, ndz)
    nlcache = nlsolver.cache.cache
    # A no-init cache is `solve!`-driven, so it has no trial-step state to reject (and no
    # `force_stop`/`nsteps` for `not_terminated` to read): a zero `ndz` there means a
    # complete inner solve returned the iterate unchanged, which is genuine convergence.
    nlcache isa NonlinearSolveNoInitCache && return false
    return (iszero(ndz) && NonlinearSolveBase.not_terminated(nlcache)) ||
        nlsolver.cache.stalled
end

"""
    compute_step!(nlsolver, integrator[, γW]) -> residual_norm

Compute one candidate nonlinear iteration and return its scaled residual or
increment norm.

# Arguments

  - `nlsolver`: the nonlinear solver whose candidate iterate and workspace are
    updated.
  - `integrator`: the differential-equation integrator that supplies the current
    state, tolerances, right-hand side, and statistics.
  - `γW`: the current `γ * dt` scaling for Newton-type methods. Fixed-point
    methods omit this argument.

# Returns

A finite, nonnegative norm when a candidate iteration was computed. Returning a
non-finite value reports divergence to [`nlsolve!`](@ref). This function mutates
the candidate iterate and its workspace and may update the integrator's nonlinear
evaluation statistics. The shared driver calls
[`OrdinaryDiffEqCore.apply_step!`](@ref) after accepting the candidate.

Solver packages that subtype [`OrdinaryDiffEqCore.AbstractNLSolver`](@ref) extend
this function to participate in the shared [`nlsolve!`](@ref) convergence loop.
"""
function compute_step! end

"""
    nlsolve!(
        nlsolver::AbstractNLSolver, integrator, cache = nothing,
        repeat_step = false
    )

Solve

```math
dt⋅f(innertmp + γ⋅z, p, t + c⋅dt) + outertmp = z
```

where `dt` is the step size and `γ` and `c` are stage constants.

# Arguments

  - `nlsolver`: a solver returned by [`build_nlsolver`](@ref), or another
    [`OrdinaryDiffEqCore.AbstractNLSolver`](@ref) implementation of this driver
    contract.
  - `integrator`: the current differential-equation integrator.
  - `cache`: the owning implicit algorithm's cache. Newton-type solvers require
    it for `W` updates; fixed-point solvers may leave it as `nothing`.
  - `repeat_step`: whether this solve is retrying the same integrator step after
    a rejection.

# Mutation and return value

The solve updates the nonlinear iterate, convergence estimate, status, failure
counters, and solver workspace. It may update the integrator state while forming
a new `W`, and its postamble updates integrator statistics and
`force_stepfail`. The return value is the result of the solver's
`SciMLBase.postamble!` method; the solvers constructed by [`build_nlsolver`](@ref)
return the converged stage increment `z`.

# Failure behavior

Non-convergence is recorded in the solver status rather than thrown: inspect it
with [`nlsolvefail`](@ref). A stale-Jacobian failure is retried once with a fresh
Jacobian. Calling a Newton-type solver without `cache` throws `ArgumentError`,
and exceptions raised by residual, Jacobian, or linear-solver evaluations
propagate.

# Extension contract

Custom nonlinear solvers extend [`compute_step!`](@ref) and [`initial_η`](@ref),
and provide the `AbstractNLSolver` state queried by the documented
`OrdinaryDiffEqCore` nonlinear-solver hooks. Candidate acceptance and finalization
are dispatched through `OrdinaryDiffEqCore.apply_step!` and `SciMLBase.postamble!`.

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
        # Checking the type, not just `nothing`: passing something else entirely
        # in this slot went unnoticed because `update_W!` happens not to read it
        # on the out-of-place path.
        cache isa Union{
            OrdinaryDiffEqCore.OrdinaryDiffEqCache,
            OrdinaryDiffEqCore.StochasticDiffEqCache,
            Nothing,
        } ||
            throw(
            ArgumentError(
                "`nlsolve!` expects the integrator cache in its third argument, got a " *
                    "$(typeof(cache))"
            )
        )
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

        if _uninformative_step(nlsolver, ndz)
            @SciMLMessage(
                lazy"Inner nonlinear solver made no progress (iter = $(iter)); the unmoved iterate is not treated as convergence",
                integrator.opts.verbose, :newton_convergence
            )
            # No new iterate exists to judge: skip the convergence and divergence
            # bookkeeping (a near-zero `ndz` would also poison the next iteration's
            # `θ`) and let the inner solver retry with its shrunken trust region. If
            # it never moves, the loop runs out and the step fails as unconverged.
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
derives it from the tolerance. The return value must be a finite nonnegative
number compatible with the integrator tolerances.

Solver packages that define an [`OrdinaryDiffEqCore.AbstractNLSolver`](@ref)
subtype extend this function when their initial estimate differs from the default
tolerance-based rule. The function does not mutate the solver or integrator and
is consumed by [`nlsolve!`](@ref).
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

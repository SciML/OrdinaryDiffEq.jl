## initialize!

@muladd function initialize!(
        nlsolver::NLSolver{<:NLNewton, false},
        integrator::SciMLBase.DEIntegrator
    )
    (; dt) = integrator
    (; cache) = nlsolver

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    nothing
end

@muladd function initialize!(
        nlsolver::NLSolver{<:NLNewton, true},
        integrator::SciMLBase.DEIntegrator
    )
    (; u, uprev, t, dt, opts) = integrator
    (; cache) = nlsolver
    (; weight) = cache

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt
    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, u,
        opts.abstol, opts.reltol, opts.internalnorm, t
    )

    nothing
end

function initialize!(
        nlsolver::NLSolver{<:NonlinearSolveAlg, false},
        integrator::SciMLBase.DEIntegrator
    )
    (; uprev, t, p, dt, opts, f) = integrator
    (; z, tmp, ztmp, γ, α, iter, cache, method, alg) = nlsolver
    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    (; ustep, tstep, k, invγdt) = cache
    # A no-init cache has no `stats` (SimpleNonlinearSolve builds every solution with
    # `stats === nothing`), and its `reinit!` only remakes the problem without evaluating
    # the residual, so neither the copies below nor the `+1` apply: nothing is counted and
    # nothing is invented. `njacs`/`nw` for a reused `W` are still recorded by
    # `_update_nlsolvealg_W_oop!` — that assembly is the integrator's own work.
    if SciMLBase.has_stats(integrator) && !(cache.cache isa NonlinearSolveNoInitCache)
        # The `reinit!` below evaluates the residual at the new `z` and *then* zeroes the
        # inner cache's counters, so that evaluation is never visible in
        # `cache.cache.stats.nf`. Left uncounted it loses one `f` call per stage.
        integrator.stats.nf += cache.cache.stats.nf + 1
        # Under `W` reuse the inner solver's "Jacobian" is the `W_ref` closure, which hands
        # back the `W` assembled here rather than evaluating anything, so its `njacs` counts
        # reads, not Jacobian evaluations; `calc_J` inside `_update_nlsolvealg_W_oop!`
        # counts the real ones. Without reuse the inner solver owns the Jacobian and its
        # count stands.
        if cache.W === nothing
            integrator.stats.njacs += cache.cache.stats.njacs
        end
        integrator.stats.nsolve += cache.cache.stats.nsolve
    end

    if cache.W !== nothing
        dtgamma = method === DIRK ? γ * dt : γ * dt / α
        W_γdt = cache.W_γdt
        first_call = iszero(W_γdt)
        # No stored `J` on this path (`W` is rebuilt from a fresh Jacobian), so a
        # dt/gamma drift past the cutoff costs a Jacobian evaluation, not just a
        # reassembly as in the in-place `new_jac`/`new_w` split.
        should_update = first_call || alg.always_new ||
            nlsolver.status === Divergence ||
            abs(inv(dtgamma) / inv(W_γdt) - 1) > oftype(dtgamma, alg.new_W_dt_cutoff)
        if should_update
            _update_nlsolvealg_W_oop!(cache, integrator, dtgamma)
            cache.new_W = true
        else
            cache.new_W = false
        end
    end

    # A no-init cache is not stepped, so there is no inner state to prime, and its
    # `reinit!` is pure heap traffic (a `remake` of the stored problem plus a re-merge of
    # the `Any`-typed kwargs, ~1.7KB per stage) spent on fields `compute_step!` no longer
    # reads: the stage problem is rebuilt there from live integrator state as a plain
    # immutable that never escapes the inner solve, so it costs no allocation. The cache
    # survives untouched as the carrier of `alg` and of `prob.f`, whose `jac` closure is
    # the reused-`W` read.
    if cache.cache isa NonlinearSolveNoInitCache
        return nothing
    end

    if f isa DAEFunction
        nlp_params = (tmp, α, tstep, invγdt, p, dt, uprev, f)
    else
        nlp_params = (tmp, γ, α, tstep, invγdt, method, p, dt, f)
    end

    refresh_stage_conditioners!(cache, tmp, γ, method, p)
    z = apply_stage_predictor!!(nlsolver)
    SciMLBase.reinit!(cache.cache, z, p = nlp_params)
    return nothing
end

function initialize!(
        nlsolver::NLSolver{<:NonlinearSolveAlg, true},
        integrator::SciMLBase.DEIntegrator
    )
    (; u, uprev, t, p, dt, opts, f) = integrator
    (; z, tmp, ztmp, γ, α, iter, cache, method, alg) = nlsolver

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    if cache.weight !== nothing
        weight = cache.weight
        calculate_residuals!(
            weight, fill!(weight, one(eltype(u))), uprev, u,
            opts.abstol, opts.reltol, opts.internalnorm, t
        )
    end

    (; ustep, atmp, tstep, k, invγdt) = cache

    # A no-init cache has no `stats` (SimpleNonlinearSolve builds every solution with
    # `stats === nothing`), and its `reinit!` only remakes the problem without evaluating
    # the residual, so neither the copies below nor the `+1` apply: nothing is counted and
    # nothing is invented. `njacs`/`nw` for the reused `W` are still recorded by
    # `_update_nlsolvealg_W!` — that assembly is the integrator's own work.
    if SciMLBase.has_stats(integrator) && !(cache.cache isa NonlinearSolveNoInitCache)
        # The `reinit!` below evaluates the residual at the new `z` and *then* zeroes the
        # inner cache's counters, so that evaluation is never visible in
        # `cache.cache.stats.nf`. Left uncounted it loses one `f` call per stage.
        integrator.stats.nf += cache.cache.stats.nf + 1
        # Under `W` reuse the inner solver's "Jacobian" is `WReuseJac`, which copies the
        # `W` assembled here rather than evaluating anything, so its `njacs` counts copies
        # (it tracks `nw`, not Jacobian evaluations). `_update_nlsolvealg_W!` counts the
        # real ones. Without reuse the inner solver owns the Jacobian and its count stands.
        if cache.W === nothing
            integrator.stats.njacs += cache.cache.stats.njacs
        end
        integrator.stats.nsolve += cache.cache.stats.nsolve
    end

    nlstep_data = f.nlstep_data
    if nlstep_data !== nothing
        atmp .= 0
        if method === COEFFICIENT_MULTISTEP
            nlstep_data.set_γ_c(nlstep_data.nlprob, (one(t), one(t), α * invγdt, tstep))
            nlstep_data.set_inner_tmp(nlstep_data.nlprob, atmp)
            nlstep_data.set_outer_tmp(nlstep_data.nlprob, tmp)
        else
            nlstep_data.set_γ_c(nlstep_data.nlprob, (dt, γ, one(t), tstep))
            nlstep_data.set_inner_tmp(nlstep_data.nlprob, tmp)
            nlstep_data.set_outer_tmp(nlstep_data.nlprob, atmp)
        end
        nlstep_data.nlprob.u0 .= @view z[nlstep_data.u0perm]
        SciMLBase.reinit!(cache.cache, nlstep_data.nlprob.u0, p = nlstep_data.nlprob.p)
    else
        if cache.W !== nothing
            dtgamma = method === DIRK ? γ * dt : γ * dt / α
            W_γdt = cache.W_γdt
            first_call = iszero(W_γdt)
            # Mirror the legacy do_newJW split: a fresh Jacobian is only needed when the
            # iteration failed or on first use, while a dt/gamma change merely requires
            # reassembling W = J - M/γdt from the stored J (jacobian2W! is O(nnz), a
            # Jacobian evaluation is not). Without adaptivity there is no error test to
            # catch the step a stale Jacobian degrades and no way to pull it back, so
            # `do_newJW` refuses to freeze `J` there and this must too.
            new_jac = first_call || alg.always_new || nlsolver.status === Divergence ||
                !integrator.opts.adaptive
            # `oftype`: `new_W_dt_cutoff` defaults to a `Rational`, and comparing a `Float64`
            # against one goes through the slow mixed-type path on every stage.
            new_w = new_jac ||
                abs(inv(dtgamma) / inv(W_γdt) - 1) > oftype(dtgamma, alg.new_W_dt_cutoff)
            if new_w
                _update_nlsolvealg_W!(cache, integrator, dtgamma, tstep, new_jac)
                cache.new_W = true
            else
                cache.new_W = false
            end
        end
        if f isa DAEFunction
            nlp_params = (tmp, ztmp, ustep, γ, α, tstep, k, invγdt, p, dt, f)
        else
            nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f)
        end
        refresh_stage_conditioners!(cache, tmp, γ, method, p)
        z = apply_stage_predictor!!(nlsolver)
        if length(get_u(cache.cache)) != length(z)
            new_prob = if cache.W !== nothing
                # W-reuse: re-point the inner jacobian at the resized W via the same mapping
                # used at build time. Only array sizes change, so the jac's concrete type —
                # and hence `cache.prob`'s concrete type — is preserved and `init` stays
                # type-stable.
                new_f = SciMLBase.remake(cache.prob.f; reuse_jac_kwargs(cache.W)...)
                SciMLBase.remake(cache.prob; f = new_f, u0 = copy(z), p = nlp_params)
            else
                SciMLBase.remake(cache.prob; u0 = copy(z), p = nlp_params)
            end
            cache.prob = new_prob
            if cache.cache isa NonlinearSolveNoInitCache
                # A no-init cache is `solve!`-driven, so it must keep terminating on its own
                # (nonzero, default) tolerances here exactly as on the build path: rebuilding
                # it with the zeroed tolerances below would leave every complete inner solve
                # returning `MaxIters`.
                cache.cache = init(
                    new_prob, cache.cache.alg;
                    verbose = integrator.opts.verbose.nonlinear_verbosity,
                    conditioning_kwargs(cache)...
                )
            else
                # Same tolerances and termination mode as `build_nlsolver`: the integrator owns
                # convergence, so the inner criterion must stay inert here too (otherwise a resized
                # cache could terminate on its own again), and the mode is part of the cache's type,
                # so rebuilding with a different one would not even be assignable.
                uT = real(eltype(z))
                cache.cache = init(
                    new_prob, cache.cache.alg;
                    verbose = integrator.opts.verbose.nonlinear_verbosity,
                    abstol = zero(uT), reltol = zero(uT),
                    termination_condition = _inner_termination(),
                    conditioning_kwargs(cache)...
                )
            end
            # re-point the estimator linsolve alias at the rebuilt inner cache
            cache.linsolve = cache.W !== nothing ? get_linear_cache(cache.cache) : nothing
        else
            # A no-init cache's `reinit!` is `remake(prob; u0, p)`: it stores the array it is
            # handed and `get_u` reads it straight back. Handing it `z` itself would alias
            # them, so after `resize!` grows `z` in place the length-mismatch rebuild above
            # could never fire and the inner problem would keep a stale `WReuseJac` over the
            # old, smaller `W` — a `SingularException` here, a silently wrong Jacobian at
            # other sizes. A no-init cache gets a copy of its own.
            znew = cache.cache isa NonlinearSolveNoInitCache ? copy(z) : z
            SciMLBase.reinit!(cache.cache, znew, p = nlp_params)
        end
    end
    return nothing
end

function _update_nlsolvealg_W!(nlcache, integrator, dtgamma, tstep, new_jac = true)
    (; J, W, uf, jac_config, du1) = nlcache
    (; f, p, uprev, alg) = integrator
    mass_matrix = f.mass_matrix
    if W isa AbstractSciMLOperator
        # Matrix-free reused W: refresh its state and gamma in place; its own `mul!`
        # supplies the Jacobian action to the inner (Krylov) solve, so there is no
        # concrete J to reassemble.
        update_coefficients!(W, uprev, p, tstep; gamma = dtgamma)
    else
        if new_jac
            if SciMLBase.has_jac(f)
                f.jac(J, uprev, p, tstep)
            elseif uf !== nothing
                uf.f = nlsolve_f(f, alg)
                uf.t = tstep
                if !(p isa SciMLBase.NullParameters)
                    uf.p = p
                end
                jacobian!(J, uf, uprev, du1, integrator, jac_config)
            end
            # `calc_J!` is bypassed here, so this is the only place the reused-`W`
            # path can record a Jacobian evaluation.
            integrator.stats.njacs += 1
        end
        jacobian2W!(W, mass_matrix, dtgamma, J)
    end
    # No estimator refresh needed: the smoothed estimate reuses the inner solver's own W
    # factorization, which the inner Newton re-factorizes itself when it refreshes W.
    nlcache.W_γdt = dtgamma
    integrator.stats.nw += 1
    return nothing
end

"""
    rescale_stale_W_rhs!(nlcache, nlsolver, integrator, nlstep_data)

Apply the stale-`W` correction of `compute_step!(::NLSolver{<:NLNewton, true}, ...)` to a
`NonlinearSolveAlg` step: when the reused `W` was assembled at a different `γΔt` than the
current step needs, the Newton direction it produces is off by that ratio.

`NLNewton` scales the increment it just computed. `NonlinearSolveAlg` cannot: the inner
solver computes *and* applies the increment inside `step!`. Scaling the right-hand side
beforehand is equivalent for the Newton descent `δu = -W \\ fu`, and unlike rewriting
`nlcache.u` afterwards it leaves the inner cache's `u`/`fu` pair mutually consistent — the
inner solver re-evaluates `fu` at the end of `step!`, so the scaling does not persist.
"""
function rescale_stale_W_rhs!(nlcache, nlsolver, integrator, nlstep_data)
    cache = nlsolver.cache
    # `W` is only reused on the plain path; otherwise the inner solver refreshes its own
    # Jacobian at every stage and there is nothing stale to correct.
    (cache.W === nothing || nlstep_data !== nothing) && return nothing
    W_γdt = cache.W_γdt
    iszero(W_γdt) && return nothing
    isdae = nlsolve_f(integrator) isa DAEFunction
    γdt = isdae ? nlsolver.α * cache.invγdt : nlsolver.γ * integrator.dt
    W_γdt ≈ γdt && return nothing
    rmul!(get_fu(nlcache), 2 / (1 + γdt / W_γdt))
    return nothing
end

"""
    inner_solve_failed(nlcache) -> Bool

Whether the inner `NonlinearSolve` cache has stopped for a reason other than reaching its
root. Stopping *successfully* is the normal path — the inner solve converged before the
integrator's own weighted test was satisfied — while every other terminal state (a collapsed
line search, an exhausted trust region, a stalled iteration) is a failure that must reach the
step controller the way `NLNewton` reports a failed linear solve. A stopped cache also makes
`step!` a no-op, so left unreported it would leave `ndz == 0` for the outer test to read as a
perfect solve.
"""
function inner_solve_failed(nlcache)
    return !NonlinearSolveBase.not_terminated(nlcache) &&
        !SciMLBase.successful_retcode(nlcache.retcode)
end

"""
    stage_unsolved(nlcache, γΔt) -> Bool

Whether the residual still says the stage equation is unsolved. `γΔt · fu` carries the units
of the iterate `u` it corrects, so `‖γΔt·fu‖ / max(‖u‖, ‖u + γΔt·fu‖)` is the relative residual
of the stage equation — invariant to the problem's units, and independent of the inner solver's
own tolerances, which `build_nlsolver` deliberately zeroes. At the cancellation floor of that
subtraction there is nothing left to solve.
"""
function stage_unsolved(nlcache, γΔt)
    fu = get_fu(nlcache)
    u = get_u(nlcache)
    resid = abs(γΔt) * maxabs(fu)
    scale = max(maxabs_axpy(γΔt, u, fu), maxabs(u))
    return resid > roundoff_level(typeof(resid)) * scale
end

# A step that moves the iterate by no more than a handful of ulps, and a residual that small
# relative to the terms it is the difference of, are both indistinguishable from zero in the
# arithmetic that produced them.
@inline roundoff_level(::Type{T}) where {T} = 8 * eps(real(one(T)))

maxabs(x::Number) = abs(x)
maxabs(x) = maximum(abs, x)

function maxabs_axpy(a, x, y)
    m = abs(a * zero(eltype(y)) + zero(eltype(x)))
    for (xi, yi) in zip(x, y)
        m = max(m, abs(xi + a * yi))
    end
    return m
end
maxabs_axpy(a, x::Number, y::Number) = abs(x + a * y)

"""
    stalled_inner_step(nlcache, ndisp, nz, fnorm_prev, γΔt) -> Bool

Whether the `step!` just taken produced no usable iterate: it moved the iterate by less than
floating-point roundoff *and* failed to reduce the residual, while the residual itself says the
stage equation is not yet solved.

The outer convergence test reads a small displacement as proximity to the root. That inference
is only sound for a step that solves the Newton system, where the displacement is `W⁻¹F` and so
vanishes only with the residual. A globalized inner solver whose line search collapses instead
returns a displacement orders of magnitude below roundoff with the residual bit-for-bit
unchanged: `ndz` is then tiny for a reason that has nothing to do with the distance to the root.

The residual has the last word. `γΔt · fu` carries the units of the iterate `u` it corrects, so
`‖γΔt·fu‖ / max(‖u‖, ‖u + γΔt·fu‖)` is the relative residual of the stage equation — invariant
to the problem's units, and independent of the inner solver's own tolerances, which
`build_nlsolver` deliberately zeroes. Only when it sits at the cancellation floor of that
subtraction is a motionless iterate genuine convergence rather than a stalled step.

`ndisp` and `nz` are `‖u - u_before‖_∞` and `‖u_before‖_∞`, and `fnorm_prev` is `‖fu‖_∞` as it
stood before the step.
"""
function stalled_inner_step(nlcache, ndisp, nz, fnorm_prev, γΔt)
    ndisp > roundoff_level(typeof(ndisp)) * nz && return false
    maxabs(get_fu(nlcache)) < fnorm_prev && return false
    return stage_unsolved(nlcache, γΔt)
end

# Scale that carries the inner residual into the units of `z`: the stage equation is
# `F(z) = 0` with `∂F/∂z = -W`, and `W` is assembled at this `γΔt` (see `initialize!`).
function residual_to_z_scale(nlsolver, isdae)
    (; α, cache, method) = nlsolver
    return (isdae || method === COEFFICIENT_MULTISTEP) ? inv(α * cache.invγdt) :
        inv(cache.invγdt)
end

function _update_nlsolvealg_W_oop!(nlcache, integrator, dtgamma)
    # Same construction as the `StaticWOperator` branch of `calc_W`: `calc_J` picks
    # user `jac` vs differentiation through `nlcache.uf` and records the evaluation
    # in `integrator.stats.njacs` — the only njacs the reused-`W` path can count.
    J_new = calc_J(integrator, nlcache)
    nlcache.W[] = J_new - integrator.f.mass_matrix * inv(dtgamma)
    nlcache.W_γdt = dtgamma
    integrator.stats.nw += 1
    return nothing
end

## compute_step!

@muladd function compute_step!(nlsolver::NLSolver{<:NonlinearSolveAlg, false}, integrator)
    (; uprev, t, p, dt, opts) = integrator
    (; z, tmp, ztmp, γ, α, cache, method) = nlsolver
    (; tstep, invγdt) = cache

    nlcache = nlsolver.cache.cache
    if nlcache isa NonlinearSolveNoInitCache
        # A no-init cache holds no iteration state, so it cannot be driven one `step!` at a
        # time: each outer iteration costs one complete inner solve, and the returned
        # solution is the only trustworthy record of it (`get_u` on the cache reads a stale
        # iterate, and `.stats`/`get_fu` do not exist). The solve bypasses the cache:
        # `solve!` on it splats `Any`-typed kwargs, a ~1KB dynamic call whose untyped
        # return poisons the rest of the branch with boxing, while a freshly built
        # problem (an immutable over the current iterate and parameter tuple) never
        # escapes the call and dispatches statically, so the solve runs allocation-free.
        # Solving from `z` also warm-starts iterations after the first at the latest outer
        # iterate. An unsuccessful solve is a genuine failure and reaches the step
        # controller the way `NLNewton` reports a failed linear solve. A complete solve
        # leaves no half-taken step behind, so there is nothing for the residual to veto.
        # The inner residual and linear-solve work is uncounted (no `stats`);
        # `_update_nlsolvealg_W_oop!` counts the Jacobian and `W` assemblies, which are
        # the integrator's own share of the work.
        if integrator.f isa DAEFunction
            nlp_params = (tmp, α, tstep, invγdt, p, dt, uprev, integrator.f)
        else
            nlp_params = (tmp, γ, α, tstep, invγdt, method, p, dt, integrator.f)
        end
        innersol = solve(NonlinearProblem(nlcache.prob.f, z, nlp_params), nlcache.alg)
        if !SciMLBase.successful_retcode(innersol.retcode)
            return convert(eltype(z), Inf)
        end
        active_u = innersol.u
        cache.stalled = false
    else
        recompute_jacobian = nlsolver.iter == 1 && (cache.W === nothing || cache.new_W)
        inner_solve_failed(nlcache) && return convert(eltype(z), Inf)
        fnorm_prev = maxabs(get_fu(nlcache))
        γΔt = residual_to_z_scale(nlsolver, nlsolve_f(integrator) isa DAEFunction)
        step!(nlcache; recompute_jacobian)
        # `step!` can *land* in a terminal state as well as start in one. Checking only on
        # entry defers a failed step to the next call, which never comes: the iterate it
        # leaves behind is unmoved, so the outer displacement test accepts the stage first.
        # An inner solver can also give up *because* there is nothing left to correct — a
        # line search collapses on a stage that is already solved — so the residual, not the
        # retcode alone, has the say.
        inner_solve_failed(nlcache) && stage_unsolved(nlcache, γΔt) &&
            return convert(eltype(z), Inf)
        active_u = get_u(nlcache)
        cache.stalled = stalled_inner_step(
            nlcache, maxabs(z .- active_u), maxabs(z), fnorm_prev, γΔt
        )
    end
    nlsolver.ztmp = active_u

    ustep = compute_ustep(tmp, γ, z, method)
    atmp = calculate_residuals(
        z .- active_u, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t
    )
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end
    return ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NonlinearSolveAlg, true}, integrator)
    (; uprev, t, p, dt, opts) = integrator
    (; z, tmp, ztmp, γ, α, cache, method) = nlsolver
    (; tstep, invγdt, atmp, ustep) = cache

    nlstep_data = integrator.f.nlstep_data
    nlcache = nlsolver.cache.cache
    if nlcache isa NonlinearSolveNoInitCache
        # A no-init cache holds no iteration state, so it cannot be driven one `step!` at a
        # time: `solve!` runs the complete inner solve and its returned solution is the only
        # trustworthy record of it (`get_u` on the cache still reads the *initial* iterate,
        # and `.stats`/`get_fu` do not exist). Each outer iteration therefore costs one full
        # inner solve; an unsuccessful one is a genuine failure and reaches the step
        # controller the way `NLNewton` reports a failed linear solve. No stale-`W` rescale
        # either: the reused `W` only serves as the inner Newton's Jacobian (`WReuseJac`),
        # where staleness costs convergence rate, not the root the full solve lands on. A
        # complete solve leaves no half-taken step behind, so nothing for the residual to veto.
        innersol = solve!(nlcache)
        if !SciMLBase.successful_retcode(innersol.retcode)
            return convert(eltype(atmp), Inf)
        end
        cache.stalled = false
    else
        recompute_jacobian = nlsolver.iter == 1 && (cache.W === nothing || cache.new_W)
        inner_solve_failed(nlcache) && return convert(eltype(atmp), Inf)
        # Read before `rescale_stale_W_rhs!`, which rewrites `fu` in place: the unscaled
        # residual is the one `step!` will produce again, so it is the like-for-like
        # comparison.
        fnorm_prev = maxabs(get_fu(nlcache))
        γΔt = residual_to_z_scale(nlsolver, nlsolve_f(integrator) isa DAEFunction)
        rescale_stale_W_rhs!(nlcache, nlsolver, integrator, nlstep_data)
        step!(nlcache; recompute_jacobian)
        # `step!` can *land* in a terminal state as well as start in one. Checking only on
        # entry defers a failed step to the next call, which never comes: the iterate it
        # leaves behind is unmoved, so the outer displacement test accepts the stage first.
        # An inner solver can also give up *because* there is nothing left to correct — a
        # line search collapses on a stage that is already solved — so the residual, not the
        # retcode alone, has the say.
        inner_solve_failed(nlcache) && stage_unsolved(nlcache, γΔt) &&
            return convert(eltype(atmp), Inf)
    end

    if nlstep_data !== nothing
        if nlcache isa NonlinearSolveNoInitCache
            # `f.nlstep_data !== nothing` selects `f.nlstep_data.nlprob`, which is `init`ed
            # with whatever inner algorithm was given, so a SimpleNonlinearSolve one lands
            # in a no-init cache here too. The solution `solve!` already returned is used
            # as-is instead of rebuilding one from `nlcache.retcode`/`nlcache.stats`/
            # `get_fu`, none of which exist on this cache.
            active_fu = innersol.resid
            nlstepsol = innersol
        else
            active_fu = get_fu(nlcache)
            # No `trace` is attached: this solution only feeds `nlprobmap` right below, and
            # polyalgorithm caches keep the trace on the active branch with no accessor for it.
            nlstepsol = SciMLBase.build_solution(
                nlcache.prob, nlcache.alg, get_u(nlcache), active_fu;
                nlcache.retcode, nlcache.stats
            )
        end
        nlstep_data.nlprobmap(ztmp, nlstepsol)
        ustep = compute_ustep!(ustep, tmp, γ, z, method)
        atmp_sub = @view(atmp[nlstep_data.u0perm])
        calculate_residuals!(
            atmp_sub, active_fu,
            @view(uprev[nlstep_data.u0perm]),
            @view(ustep[nlstep_data.u0perm]), opts.abstol,
            opts.reltol, opts.internalnorm, t
        )
        # Norm must be computed over the nlstep subset only.
        # The full-length `atmp` has zeros in non-`u0perm` slots (set in
        # `initialize!`), which would otherwise dilute the scalar norm
        # (e.g. `ODE_DEFAULT_NORM` divides by `length`), causing the Newton
        # convergence check to declare success ~sqrt(n_full/n_sub) early.
        ndz = opts.internalnorm(atmp_sub, t)
        # This branch already measures the residual rather than a displacement, so there is
        # no displacement to be misled by.
        cache.stalled = false
    else
        active_u = nlcache isa NonlinearSolveNoInitCache ? innersol.u : get_u(nlcache)
        @.. broadcast = false ztmp = active_u
        ustep = compute_ustep!(ustep, tmp, γ, z, method)
        @.. broadcast = false atmp = z - ztmp
        if !(nlcache isa NonlinearSolveNoInitCache)
            cache.stalled = stalled_inner_step(
                nlcache, maxabs(atmp), maxabs(z), fnorm_prev, γΔt
            )
        end
        calculate_residuals!(
            atmp, atmp, uprev, ustep, opts.abstol, opts.reltol,
            opts.internalnorm, t
        )
        ndz = opts.internalnorm(atmp, t)
    end

    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end
    ndz
end

"""
    compute_step!(nlsolver::NLSolver{<:NLNewton}, integrator)

Compute next iterate of numerically stable modified Newton iteration
that is specialized for implicit methods.

Please check
https://github.com/SciML/DiffEqDevMaterials/blob/master/newton/output/main.pdf
for more details.

# References

M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8).

Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).
"""
@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, false}, integrator, γW)
    (; uprev, t, p, dt, opts) = integrator
    (; z, tmp, ztmp, γ, α, cache, method) = nlsolver
    (; tstep, W, invγdt) = cache

    f = nlsolve_f(integrator)

    if f isa DAEFunction
        _uprev = get_dae_uprev(integrator, uprev)
        ztmp, ustep = _compute_rhs(tmp, α, tstep, invγdt, p, _uprev, f, z)
    else
        ztmp, ustep = _compute_rhs(tmp, γ, α, tstep, invγdt, method, p, dt, f, z)
    end

    if SciMLBase.has_stats(integrator)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    # update W
    if W isa Union{WOperator, StaticWOperator}
        update_coefficients!(W, ustep, p, tstep)
    elseif W isa AbstractSciMLOperator
        error(
            "Out-of-place NLNewton does not support non-concrete W of type $(typeof(W)). " *
                "Use the in-place form (define f!(du, u, p, t)) or supply a concrete jac_prototype."
        )
    end
    dz = _restructure_state(ztmp, W \ _vec(ztmp))
    dz = relax(dz, nlsolver, integrator, f)
    if SciMLBase.has_stats(integrator)
        integrator.stats.nsolve += 1
    end

    atmp = calculate_residuals(
        dz, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t
    )
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    nlsolver.ztmp = z .- dz

    ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, true}, integrator, γW)
    (; uprev, t, p, dt, opts) = integrator
    (; z, tmp, ztmp, γ, α, iter, cache, method) = nlsolver
    (; W_γdt, ustep, tstep, k, atmp, dz, W, new_W, invγdt, linsolve, weight) = cache

    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    if SciMLBase.has_stats(integrator)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    if isdae
        _uprev = get_dae_uprev(integrator, uprev)
        b, ustep = _compute_rhs!(tmp, ztmp, ustep, α, tstep, k, invγdt, p, _uprev, f, z)
    else
        b, ustep = _compute_rhs!(
            tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f, z
        )
    end

    # update W
    if W isa Union{WOperator, StaticWOperator}
        update_coefficients!(W, ustep, p, tstep)
    elseif W isa AbstractSciMLOperator
        # logic for generic AbstractSciMLOperator does not yet support partial state updates, so provide full state
        update_coefficients!(W, ustep, p, tstep; gamma = γW, transform = true)
    end

    if integrator.opts.adaptive
        reltol = integrator.opts.reltol
    else
        reltol = eps(real(one(eltype(dz))))
    end

    if is_always_new(nlsolver) || (iter == 1 && new_W)
        linres = dolinsolve(
            integrator, linsolve; A = W, b = _vec(b), linu = _vec(dz),
            reltol = reltol
        )
    else
        linres = dolinsolve(
            integrator, linsolve; A = nothing, b = _vec(b), linu = _vec(dz),
            reltol = reltol
        )
    end

    if !SciMLBase.successful_retcode(linres.retcode) &&
            linres.retcode != SciMLBase.ReturnCode.Default
        return convert(eltype(atmp), Inf)
    end

    if SciMLBase.has_stats(integrator)
        integrator.stats.nsolve += 1
    end

    # relaxed Newton
    # Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
    # Equations. A Review, by Christopher A. Kennedy and Mark H. Carpenter
    # page 54.
    γdt = isdae ? α * invγdt : γ * dt

    !(W_γdt ≈ γdt) && (rmul!(dz, 2 / (1 + γdt / W_γdt)))
    relax!(dz, nlsolver, integrator, f)

    calculate_residuals!(
        atmp, dz, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t
    )
    ndz = opts.internalnorm(atmp, t)

    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    @.. broadcast = false ztmp = z - dz

    ndz
end

function get_dae_uprev(integrator, uprev)
    # not all predictors are uprev, for other forms of predictors, defined in u₀
    return if isdefined(integrator.cache, :u₀)
        integrator.cache.u₀
    else
        uprev
    end
end

function _compute_rhs(tmp, α, tstep, invγdt, p, uprev, f::TF, z) where {TF <: DAEFunction}
    ustep = @.. uprev + z
    dustep = @. (tmp + α * z) * invγdt
    ztmp = f(dustep, ustep, p, tstep)
    return ztmp, ustep
end

function compute_ustep(tmp, γ, z, method)
    return if method === COEFFICIENT_MULTISTEP
        z
    else
        @. tmp + γ * z
    end
end

function compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        ustep = z
    else
        @.. ustep = tmp + γ * z
    end
    return ustep
end

function _compute_rhs(tmp, γ, α, tstep, invγdt, method::MethodType, p, dt, f::F, z) where {F}
    mass_matrix = f.mass_matrix
    ustep = compute_ustep(tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        # tmp = outertmp ./ hγ
        if mass_matrix === I
            ztmp = tmp .+ f(z, p, tstep) .- (α * invγdt) .* z
        else
            update_coefficients!(mass_matrix, z, p, tstep)
            ztmp = tmp .+ f(z, p, tstep) .- (mass_matrix * z) .* (α * invγdt)
        end
    else
        if mass_matrix === I
            ztmp = (dt * f(ustep, p, tstep) - z) * invγdt
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            ztmp = (dt * f(ustep, p, tstep) - mass_matrix * z) * invγdt
        end
    end
    return ztmp, ustep
end

function _compute_rhs!(
        tmp, ztmp, ustep, α, tstep, k,
        invγdt, p, uprev, f::TF, z
    ) where {TF <: DAEFunction}
    @.. broadcast = false ztmp = (tmp + α * z) * invγdt
    @.. ustep = uprev + z
    f(k, ztmp, ustep, p, tstep)
    return _vec(k), ustep
end

function _compute_rhs!(
        tmp, ztmp, ustep, γ, α, tstep, k,
        invγdt, method::MethodType, p, dt, f, z
    )
    mass_matrix = f.mass_matrix
    ustep = compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        f(k, z, p, tstep)
        if mass_matrix === I
            @.. broadcast = false ztmp = tmp + k - (α * invγdt) * z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @.. broadcast = false ztmp = tmp + k - (α * invγdt) * ztmp
        end
    else
        f(k, ustep, p, tstep)
        if mass_matrix === I
            @.. ztmp = (dt * k - z) * invγdt
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @.. ztmp = (dt * k - ztmp) * invγdt
        end
    end
    return _vec(ztmp), ustep
end

function _compute_rhs!(
        tmp::Array, ztmp::Array, ustep::Array, α, tstep, k,
        invγdt, p, uprev, f::TF, z
    ) where {TF <: DAEFunction}
    @inbounds @simd ivdep for i in eachindex(z)
        ztmp[i] = (tmp[i] + α * z[i]) * invγdt
    end
    @inbounds @simd ivdep for i in eachindex(z)
        ustep[i] = uprev[i] + z[i]
    end
    f(k, ztmp, ustep, p, tstep)

    return _vec(k), ustep
end

function _compute_rhs!(
        tmp::Array, ztmp::Array, ustep::Array, γ, α, tstep, k,
        invγdt, method::MethodType, p, dt, f, z
    )
    mass_matrix = f.mass_matrix
    ustep = compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        f(k, z, p, tstep)
        if mass_matrix === I
            @inbounds @simd ivdep for i in eachindex(z)
                ztmp[i] = tmp[i] + k[i] - (α * invγdt) * z[i]
            end
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))

            @inbounds @simd ivdep for i in eachindex(z)
                ztmp[i] = tmp[i] + k[i] - (α * invγdt) * ztmp[i]
            end
        end
    else
        f(k, ustep, p, tstep)
        if mass_matrix === I
            @inbounds @simd ivdep for i in eachindex(z)
                ztmp[i] = (dt * k[i] - z[i]) * invγdt
            end
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @inbounds @simd ivdep for i in eachindex(z)
                ztmp[i] = (dt * k[i] - ztmp[i]) * invγdt
            end
        end
    end

    return _vec(ztmp), ustep
end

## relax!
function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF}
    return relax!(dz, nlsolver, integrator, f, relax(nlsolver))
end
function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF}
    return relax(dz, nlsolver, integrator, f, relax(nlsolver))
end
function relax!(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Nothing
    ) where {TF}
    return dz
end
function relax!(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Number
    ) where {TF}
    if !iszero(r)
        rmul!(dz, 1 - r)
    end
    return dz
end

function relax!(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        linesearch
    ) where {TF}
    let dz = dz,
            integrator = integrator,
            nlsolver = nlsolver,
            f = f,
            linesearch = linesearch

        (; uprev, t, p, dt, opts, isdae) = integrator
        (; z, tmp, ztmp, γ, iter, α, cache, method) = nlsolver
        (; ustep, atmp, tstep, k, invγdt) = cache
        function resid(z)
            # recompute residual (rhs)
            if isdae
                _uprev = get_dae_uprev(integrator, uprev)
                b, ustep2 = _compute_rhs!(
                    tmp, ztmp, ustep, α, tstep, k, invγdt, p, _uprev, f::TF, z
                )
            else
                b, ustep2 = _compute_rhs!(
                    tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f, z
                )
            end
            calculate_residuals!(
                atmp, b, uprev, ustep2, opts.abstol, opts.reltol,
                opts.internalnorm, t
            )
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            local z = @.. atmp = nlsolver.z - dz * α
            res = resid(z)
            return res
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α + ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            ∂ϕ∂α = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, ∂ϕ∂α
        end
        α0 = one(eltype(ustep))
        ϕ0, dϕ0 = ϕdϕ(zero(α0))
        α, _ = linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        @.. dz = dz * α
        return dz
    end
end

function relax(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Number
    ) where {TF}
    if !iszero(r)
        dz = (1 - r) * dz
    end
    return dz
end

function relax(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Nothing
    ) where {TF}
    return dz
end

function relax(
        dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        linesearch
    ) where {TF}
    let dz = dz,
            integrator = integrator,
            nlsolver = nlsolver,
            f = f,
            linesearch = linesearch

        (; uprev, t, p, dt, opts) = integrator
        (; z, tmp, ztmp, γ, iter, cache, method) = nlsolver
        (; ustep, atmp, tstep, k, invγdt) = cache
        function resid(z)
            # recompute residual (rhs)
            if f isa DAEFunction
                _uprev = get_dae_uprev(integrator, uprev)
                ztmp, ustep2 = _compute_rhs(tmp, α, tstep, invγdt, p, dt, _uprev, f, z)
            else
                ztmp, ustep2 = _compute_rhs(tmp, γ, α, tstep, invγdt, method, p, f, z)
            end
            atmp = calculate_residuals(
                b, uprev, ustep2, opts.abstol, opts.reltol,
                opts.internalnorm, t
            )
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            local z = @.. nlsolver.z - dz * α
            return resid(z)
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α + ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            ∂ϕ∂α = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, ∂ϕ∂α
        end
        α0 = one(eltype(dz))
        ϕ0, dϕ0 = ϕdϕ(zero(α0))
        α, _ = linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        dz = dz * α
        return dz
    end
end

## resize!

function Base.resize!(nlcache::NLNewtonCache, ::AbstractNLSolver, integrator, i::Int)
    resize!(nlcache.ustep, i)
    resize!(nlcache.k, i)
    resize!(nlcache.atmp, i)
    resize!(nlcache.dz, i)
    resize!(nlcache.du1, i)

    resize_jac_config!(nlcache, integrator)
    resize!(nlcache.weight, i)

    # resize J and W (or rather create new ones of appropriate size and type)
    resize_J_W!(nlcache, integrator, i)

    return nothing
end

function Base.resize!(nlcache::NonlinearSolveCache, ::AbstractNLSolver, integrator, i::Int)
    nlcache.ustep === nothing || resize!(nlcache.ustep, i)
    nlcache.k === nothing || resize!(nlcache.k, i)
    nlcache.atmp === nothing || resize!(nlcache.atmp, i)
    nlcache.du1 === nothing || resize!(nlcache.du1, i)
    nlcache.weight === nothing || resize!(nlcache.weight, i)
    nlcache.dz === nothing || resize!(nlcache.dz, i)
    nlcache.jac_config === nothing || resize_jac_config!(nlcache, integrator)
    # The reused out-of-place `W` is a `Ref` over a static matrix, which has no `resize!`.
    nlcache.W === nothing || nlcache.W isa Ref || resize_J_W!(nlcache, integrator, i)
    resize_conditioner!(nlcache.precondition, i)
    resize_conditioner!(nlcache.postcondition, i)
    # `nlcache.linsolve` is re-pointed by `initialize!` when it rebuilds the inner cache on the
    # next length mismatch, before any estimator solve — nothing to rebuild here.
    nlcache.W_γdt = zero(nlcache.W_γdt)
    nlcache.new_W = true
    return nothing
end

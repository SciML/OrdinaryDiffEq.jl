# ── MREIL: Jacobian and W matrix of the fast component ────────────────────────
#
# `cache.jac_f` is what `alg_cache` handed `build_J_W`: the fast component with the user
# `jac` and the `jac_prototype`/`sparsity` resolved onto it from wherever they were
# written. It is by construction the description `J` and `W` were built from, so the
# `jac_prototype` below is the merged one and cannot drift from the structure `J` has.
#
# Read the Jacobian only from there. Re-deriving it from `integrator.f` reads either the
# unmerged `f1`, whose fields are `nothing` in exactly the outer-placement configuration
# this exists to support, or the outer function, whose `jac` is paired with the outer
# `jac_prototype` rather than the merged one. Because `_mreil_jac_function` resolved
# `jac` once, `has_jac(cache.jac_f)` agrees with `has_jac(integrator.f)` and the
# fall-through below is reachable only with no user Jacobian anywhere — where `calc_J!`
# and `calc_J` take their AD branch against `nlsolve_f(f, ::MREIL) == f1`.

# Give `J` `jac_prototype`'s sparsity structure with all stored values zeroed, ready for
# a user-supplied `f.jac` to write into: `f.jac` only assigns into already-stored entries,
# so J's structure must track jac_prototype's rather than whatever `f.jac` fills in.
# https://github.com/SciML/OrdinaryDiffEq.jl/issues/2653
_mreil_prepare_jac!(J, jac_prototype) = nothing
function _mreil_prepare_jac!(J::SparseMatrixCSC, jac_prototype::SparseMatrixCSC)
    if !(
            size(J) == size(jac_prototype) && getcolptr(J) == getcolptr(jac_prototype) &&
                rowvals(J) == rowvals(jac_prototype)
        )
        # A structural copy, not a broadcast: broadcasting prunes numerical zeros, which
        # would drop entries the prototype stores as zero — and working around that by
        # writing nonzero values into `jac_prototype` first would corrupt the user's
        # array, which `remake` shares across solves.
        copyto!(J, jac_prototype)
    end
    nonzeros(J) .= false
    return nothing
end

function _mreil_calc_J!(J, integrator, cache)
    jac_f = cache.jac_f
    SciMLBase.has_jac(jac_f) || return calc_J!(J, integrator, cache)
    _mreil_prepare_jac!(J, jac_f.jac_prototype)
    jac_f.jac(J, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.njacs += 1
    return nothing
end

function _mreil_calc_J(integrator, cache)
    jac_f = cache.jac_f
    SciMLBase.has_jac(jac_f) || return calc_J(integrator, cache)
    integrator.stats.njacs += 1
    return jac_f.jac(integrator.uprev, integrator.p, integrator.t)
end

# J is taken at (uprev, t), and the only thing that may reuse it is a retry of a
# rejected attempt at the very same step — nothing else is guaranteed to leave `p`, the
# problem functions and the state alone. `repeat_step` is never passed as `true` for a
# non-composite algorithm, so the retry has to be recognised from the integrator.
#
# `success_iter` counts steps the core loop has applied and `iter` counts attempts; both
# are bumped in `loopheader!` before `perform_step!` runs, and `reinit!` zeroes both
# whether or not it re-runs `initialize!`. So `success_iter` unchanged *and* `iter` one
# past the last attempt holds exactly on a retry. Every way a user can touch a live
# integrator — a callback, `set_u!`, `set_t!`, mutating `p`, `reinit!`, or anything
# between two `step!` calls — runs only after an accepted step or after `reinit!`, and
# both break the test. Keying on `(uprev, t)` instead did not: `reinit!` restores them.
function _mreil_J_isreusable(integrator, cache)
    cache.J_isvalid || return false
    cache.J_success_iter == integrator.success_iter || return false
    return cache.J_iter + 1 == integrator.iter
end

# Refresh the concrete Jacobian behind `cache.W`, if there is one: a `WOperator` over
# a matrix-free JVP or over a linear-operator `f1` carries no matrix to fill in, and
# W itself is then re-evaluated lazily from the operator on every use.
function _mreil_update_J!(integrator, cache)
    W = cache.W
    J = W isa WOperator ? W.J : cache.J
    if J isa AbstractMatrix && !_mreil_J_isreusable(integrator, cache)
        _mreil_calc_J!(J, integrator, cache)
        cache.J_isvalid = true
    end
    # Stamp the attempt even when J was reused, so that a step rejected twice in a row
    # keeps reusing it instead of falling outside the one-attempt window on the retry.
    cache.J_iter = integrator.iter
    cache.J_success_iter = integrator.success_iter
    return nothing
end

# Build W = J - M/h_fast for one extrapolation column. Mirrors `calc_W!`'s dispatch:
# an operator W only needs its `gamma` (and its `(u, p, t)`, for a time-dependent
# operator) updated, while a concrete W is written from the concrete J. `nw` counts as
# `calc_W`'s does: every branch that produces a new W — lazy or concrete — except a
# user-supplied `W_prototype` operator, where only its coefficients move and no W is
# built at all.
function _mreil_update_W!(integrator, cache, h_fast)
    (; uprev, p, t, f) = integrator
    W = cache.W
    if W isa WOperator
        update_coefficients!(W, uprev, p, t; gamma = h_fast)
        if W.J isa AbstractMatrix
            jacobian2W!(W._concrete_form, f.mass_matrix, h_fast, W.J)
        end
        integrator.stats.nw += 1
    elseif W isa AbstractSciMLOperator
        update_coefficients!(W, uprev, p, t; gamma = h_fast)
    else
        jacobian2W!(W, f.mass_matrix, h_fast, cache.J)
        integrator.stats.nw += 1
    end
    return nothing
end

# W = J - M/h_fast, i.e. -(M - h_fast*J)/h_fast. `jacobian2W` needs an `AbstractMatrix`
# J, so a scalar state gets its own methods; the default `M = I` has `λ = true` and
# reduces to `J - inv(h_fast)`.
@inline function _mreil_W(mass_matrix::UniformScaling, h_fast, J::Number)
    return J - mass_matrix.λ * inv(h_fast)
end
@inline _mreil_W(mass_matrix::Number, h_fast, J::Number) = J - mass_matrix * inv(h_fast)
@inline function _mreil_W(mass_matrix, h_fast, J::Number)
    throw(
        ArgumentError(
            "MREIL cannot apply a mass matrix of type $(typeof(mass_matrix)) to a " *
                "scalar state; use a `UniformScaling` or a scalar mass matrix, or make " *
                "the state an array."
        )
    )
end
@inline _mreil_W(mass_matrix, h_fast, J) = jacobian2W(mass_matrix, h_fast, J)

@inline _mreil_factorize(W::Number) = W
@inline _mreil_factorize(W) = DiffEqBase.default_factorize(W)

@inline _mreil_restructure(template::Number, x) = oftype(template, x)
@inline _mreil_restructure(template, x) = _reshape(x, size(template))

# ── MREIL initialize! ─────────────────────────────────────────────────────────

function initialize!(integrator, cache::MREILCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)  # f1
    integrator.stats.nf2 += 1  # f2
    return integrator.fsalfirst .+= cache.tmp
end

function initialize!(integrator, cache::MREILConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

# ── MREIL perform_step! (in-place, MutableCache) ──────────────────────────────
#
# Same skeleton as MREEF — nj macro intervals with a frozen slow rate, m fast
# substeps each, then Aitken–Neville extrapolation over T[1..order] — but the
# fast substep is a linearly implicit Euler step
#   (M - h_fast*J) * Δ = h_fast*(k_slow + k_fast),  T[j] += Δ
# with J the Jacobian of the fast component f1, where the stiffness lives.

function perform_step!(integrator, cache::MREILCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, weight, k_slow, k_fast, linsolve_tmp, T, W) = cache
    alg = unwrap_alg(integrator, true)
    m = alg.m
    order = alg.order
    ns = cache.ns

    # One Jacobian per macro step, taken at (uprev, t) and shared by every column
    # of the extrapolation table. The frozen J is part of the base method, so a
    # column-dependent J would give each T[j] a different asymptotic h-expansion
    # and Richardson extrapolation would no longer cancel the leading error
    # terms. Refreshing J once per step is also what the single-rate extrapolated
    # linearly implicit Euler (SEULEX) does.
    _mreil_update_J!(integrator, cache)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        # h_fast depends only on the column index, so W is constant over the whole
        # column: build and factorize it once here, then every one of the nj*m
        # solves below reuses that factorization (`dolinsolve` without `A`). W is a
        # raw matrix or operator here — the factorization lives inside the LinearSolve
        # cache — so singularity only surfaces through the solve's retcode below.
        _mreil_update_W!(integrator, cache, h_fast)
        refactorize = true

        Tj = T[j]
        vecTj = _vec(Tj)
        @.. broadcast = false Tj = uprev

        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac

            # Slow evaluation (f2): frozen for all m fast substeps
            f.f2(k_slow, Tj, p, t_mac)
            integrator.stats.nf2 += 1

            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                f.f1(k_fast, Tj, p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                # first one has already evaluated the whole right-hand side there —
                # the left endpoint the Hermite interpolant needs. MREIL is not FSAL,
                # so nothing else refreshes k[1] between steps.
                if j == 1 && i_mac == 1 && i_fast == 1
                    @.. broadcast = false integrator.fsalfirst = k_slow + k_fast
                end

                @.. broadcast = false linsolve_tmp = -(k_slow + k_fast)
                linres = if refactorize
                    refactorize = false
                    dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp))
                else
                    dolinsolve(integrator, cache.linsolve; b = _vec(linsolve_tmp))
                end
                integrator.stats.nsolve += 1
                # LinearSolve's factorizations report a singular W as a `Failure`
                # retcode, but an algorithm may leave the retcode at `Default`, where
                # a failure is only visible through non-finite solution values. Either
                # way it is a property of the linear solver at this h_fast, not of the
                # step error, so fail the attempt the way `nlsolve!` does: `adaptive`
                # then retries at a smaller dt, and a fixed-dt solve aborts loudly
                # instead of accepting a `u` that was never written.
                if !(
                        SciMLBase.successful_retcode(linres.retcode) ||
                            (
                            linres.retcode == SciMLBase.ReturnCode.Default &&
                                all(isfinite, linres.u)
                        )
                    )
                    integrator.force_stepfail = true
                    return nothing
                end

                vecdz = _vec(linres.u)
                @.. broadcast = false vecTj = vecTj + vecdz
            end
        end
    end

    # Aitken–Neville Richardson extrapolation (in-place, reverse-row order). The
    # ratio is kept exact as a `Rational` so the division below rounds once in the
    # state's eltype — an Int/Int ratio would cap a BigFloat solve at Float64 accuracy.
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] // ns[j - k]
            @.. broadcast = false tmp = (T[j] - T[j - 1]) / (ratio - 1)
            @.. broadcast = false T[j] = T[j] + tmp
        end
    end

    @.. broadcast = false u = T[order]

    if integrator.opts.adaptive
        @.. broadcast = false tmp = T[order] - T[order - 1]
        calculate_residuals!(
            atmp,
            tmp,
            uprev,
            u,
            integrator.opts.abstol,
            integrator.opts.reltol,
            integrator.opts.internalnorm,
            t,
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    f.f1(integrator.fsallast, u, p, t + dt)
    f.f2(k_slow, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    @.. broadcast = false integrator.fsallast = integrator.fsallast + k_slow
    return nothing
end

# ── MREIL perform_step! (out-of-place, ConstantCache) ─────────────────────────

@muladd function perform_step!(integrator, cache::MREILConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, true)
    m = alg.m
    order = alg.order
    ns = cache.ns
    T = cache.T
    mass_matrix = f.mass_matrix

    J = _mreil_calc_J(integrator, cache)

    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        W = _mreil_factorize(_mreil_W(mass_matrix, h_fast, J))
        integrator.stats.nw += 1
        if !issuccess_W(W)
            # A singular W is a property of this h_fast, not of the step error: fail
            # the attempt the way `nlsolve!` does, so `adaptive` retries at a smaller
            # dt and a fixed-dt solve aborts loudly instead of accepting a `u` that
            # was never written.
            integrator.force_stepfail = true
            return nothing
        end

        u_cur = uprev
        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac
            k_slow = f.f2(u_cur, p, t_mac)
            integrator.stats.nf2 += 1
            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                k_fast = f.f1(u_cur, p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                # first one has already evaluated the whole right-hand side there —
                # the left endpoint the Hermite interpolant needs. MREIL is not FSAL,
                # so nothing else refreshes k[1] between steps.
                if j == 1 && i_mac == 1 && i_fast == 1
                    integrator.fsalfirst = @.. broadcast = false k_slow + k_fast
                end

                rhs = @.. broadcast = false -(k_slow + k_fast)
                dz = _mreil_restructure(uprev, W \ _vec(rhs))
                integrator.stats.nsolve += 1
                u_cur = @.. broadcast = false u_cur + dz
            end
        end
        T[j] = u_cur
    end

    # Aitken–Neville Richardson extrapolation. The ratio is kept exact as a
    # `Rational` so the division below rounds once in the state's eltype — an
    # Int/Int ratio would cap a BigFloat solve at Float64 accuracy.
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] // ns[j - k]
            T[j] = @.. broadcast = false T[j] + (T[j] - T[j - 1]) / (ratio - 1)
        end
    end

    integrator.u = T[order]

    if integrator.opts.adaptive
        utilde = @.. broadcast = false T[order] - T[order - 1]
        atmp = calculate_residuals(
            utilde,
            uprev,
            integrator.u,
            integrator.opts.abstol,
            integrator.opts.reltol,
            integrator.opts.internalnorm,
            t,
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    integrator.fsallast = f.f1(integrator.u, p, t + dt) + f.f2(integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function initialize!(integrator, cache::MRIGARKImplicitCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst .+= cache.tmp
end

function initialize!(integrator, cache::MRIGARKImplicitConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::MRIGARKImplicitCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, z, fS, zemb, nlsolver, tab) = cache
    (; Δc, W0, W1, Wemb0, Wemb1, γ0, q) = tab

    f.f1(integrator.fsalfirst, uprev, p, t)
    f.f2(cache.f1eval, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    @.. broadcast = false integrator.fsalfirst = integrator.fsalfirst + cache.f1eval

    alg = unwrap_alg(integrator, true)
    m = alg.m
    s = length(Δc)
    stats = integrator.stats

    markfirststage!(nlsolver)

    @.. broadcast = false z[1] = uprev
    cprev = zero(eltype(Δc))
    for i in 1:s
        f.f2(fS[i], z[i], p, t + cprev * dt)
        stats.nf2 += 1
        cnext = cprev + Δc[i]
        if iszero(γ0[i])
            _mrigark_substage!(
                z[i + 1], z[i], cache, f, p, t, dt, Δc[i], cprev, m, q,
                view(W0, i, :), view(W1, i, :), fS, i, stats
            )
        else
            @.. broadcast = false nlsolver.tmp = z[i]
            for j in 1:i
                ω̄ = W0[i, j] + W1[i, j] / 2
                iszero(ω̄) && continue
                @.. broadcast = false nlsolver.tmp = nlsolver.tmp + dt * ω̄ * fS[j]
            end
            nlsolver.c = cnext
            @.. broadcast = false nlsolver.z = dt * fS[i]
            w = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            @.. broadcast = false z[i + 1] = nlsolver.tmp + nlsolver.γ * w
        end
        cprev = cnext
    end
    @.. broadcast = false u = z[s + 1]

    if integrator.opts.adaptive
        if isempty(Wemb0)
            @.. broadcast = false tmp = u - z[s]
        else
            w1emb = isempty(Wemb1) ? nothing : Wemb1
            _mrigark_substage!(
                zemb, z[s], cache, f, p, t, dt, Δc[s], cprev - Δc[s], m, q,
                Wemb0, w1emb, fS, s, stats
            )
            @.. broadcast = false tmp = u - zemb
        end
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    f.f1(integrator.fsallast, u, p, t + dt)
    f.f2(cache.f1eval, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    @.. broadcast = false integrator.fsallast = integrator.fsallast + cache.f1eval
    return nothing
end

@muladd function perform_step!(
        integrator, cache::MRIGARKImplicitConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; nlsolver, tab, z, fS) = cache
    (; Δc, W0, W1, Wemb0, Wemb1, γ0, q) = tab

    integrator.fsalfirst = f.f1(uprev, p, t) + f.f2(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1

    alg = unwrap_alg(integrator, true)
    m = alg.m
    s = length(Δc)
    stats = integrator.stats

    markfirststage!(nlsolver)

    z[1] = uprev
    cprev = zero(eltype(Δc))
    for i in 1:s
        fS[i] = f.f2(z[i], p, t + cprev * dt)
        stats.nf2 += 1
        cnext = cprev + Δc[i]
        if iszero(γ0[i])
            z[i + 1] = _mrigark_substage(
                z[i], f, p, t, dt, Δc[i], cprev, m, q,
                view(W0, i, :), view(W1, i, :), fS, i, stats
            )
        else
            tmp = z[i]
            for j in 1:i
                ω̄ = W0[i, j] + W1[i, j] / 2
                iszero(ω̄) && continue
                tmp = @.. broadcast = false tmp + dt * ω̄ * fS[j]
            end
            nlsolver.tmp = tmp
            nlsolver.c = cnext
            nlsolver.z = dt * fS[i]
            w = nlsolve!(nlsolver, integrator, cache, repeat_step)
            nlsolvefail(nlsolver) && return
            z[i + 1] = @.. broadcast = false nlsolver.tmp + nlsolver.γ * w
        end
        cprev = cnext
    end
    integrator.u = z[s + 1]

    if integrator.opts.adaptive
        utilde = if isempty(Wemb0)
            @.. broadcast = false integrator.u - z[s]
        else
            w1emb = isempty(Wemb1) ? nothing : Wemb1
            zemb = _mrigark_substage(
                z[s], f, p, t, dt, Δc[s], cprev - Δc[s], m, q,
                Wemb0, w1emb, fS, s, stats
            )
            @.. broadcast = false integrator.u - zemb
        end
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    integrator.fsallast = f.f1(integrator.u, p, t + dt) + f.f2(integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

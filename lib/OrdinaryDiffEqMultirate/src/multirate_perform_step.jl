# ── Extrapolated multirate methods: step-count sequence ───────────────────────

@inline function _extrapolation_sequence(seq::Symbol, order::Int)
    if seq === :harmonic
        return ntuple(j -> j, order)
    elseif seq === :romberg
        return ntuple(j -> 1 << (j - 1), order)
    else
        throw(ArgumentError("unknown sequence `$seq`, choose :harmonic or :romberg"))
    end
end

# ── MREEF initialize! ─────────────────────────────────────────────────────────

function initialize!(integrator, cache::MREEFCache)
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

function initialize!(integrator, cache::MREEFConstantCache)
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

# ── MREEF perform_step! (in-place, MutableCache) ──────────────────────────────
#
# Base multirate Euler with nj macro intervals, m fast substeps each:
#   1. k_slow = f.f2(u, p, t_mac)  — frozen slow rate for the macro interval
#   2. m fast substeps: u += h_fast*(k_slow + f.f1(u, p, t_fast))
# f1 = fast/stiff (large eigenvalues), f2 = slow/non-stiff (SciML convention).
# Then apply Aitken–Neville Richardson extrapolation over T[1..order].

function perform_step!(integrator, cache::MREEFCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, k_slow, k_fast, T) = cache
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _extrapolation_sequence(alg.seq, order)

    # Fill first tableau column: T[j] = base method with ns[j] macro intervals
    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        @.. broadcast = false T[j] = uprev

        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac

            # Slow evaluation (f2): frozen for all m fast substeps
            f.f2(k_slow, T[j], p, t_mac)
            integrator.stats.nf2 += 1

            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                f.f1(k_fast, T[j], p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                @.. broadcast = false T[j] = T[j] + h_fast * k_slow + h_fast * k_fast
            end
        end
    end

    # Aitken–Neville Richardson extrapolation (in-place, reverse-row order)
    # Formula: T[j] <- T[j] + (T[j] - T[j-1]) / (ns[j]/ns[j-k] - 1)
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
            @.. broadcast = false tmp = (T[j] - T[j - 1]) / (ratio - 1)
            @.. broadcast = false T[j] = T[j] + tmp
        end
    end

    @.. broadcast = false u = T[order]

    return if integrator.opts.adaptive
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
end

# ── MREEF perform_step! (out-of-place, ConstantCache) ─────────────────────────

@muladd function perform_step!(integrator, cache::MREEFConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    m = alg.m
    order = alg.order
    ns = _extrapolation_sequence(alg.seq, order)
    T = cache.T

    for j in 1:order
        nj = ns[j]
        h_mac = dt / nj
        h_fast = h_mac / m

        u_cur = uprev
        for i_mac in 1:nj
            t_mac = t + (i_mac - 1) * h_mac
            k_slow = f.f2(u_cur, p, t_mac)
            integrator.stats.nf2 += 1
            for i_fast in 1:m
                t_fast = t_mac + (i_fast - 1) * h_fast
                k_fast = f.f1(u_cur, p, t_fast)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                u_cur = @.. broadcast = false u_cur + h_fast * k_slow + h_fast * k_fast
            end
        end
        T[j] = u_cur
    end

    # Aitken–Neville Richardson extrapolation
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
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
end

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
        # `jac_prototype`'s stored values must be made nonzero first: the broadcast below
        # prunes numerical zeros, which would drop those entries from `J`.
        nonzeros(jac_prototype) .= true
        J .= true .* jac_prototype
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
# operator) updated, while a concrete W is written from the concrete J.
function _mreil_update_W!(integrator, cache, h_fast)
    (; uprev, p, t, f) = integrator
    W = cache.W
    if W isa WOperator
        update_coefficients!(W, uprev, p, t; gamma = h_fast)
        if W.J isa AbstractMatrix
            jacobian2W!(W._concrete_form, f.mass_matrix, h_fast, W.J)
        end
    elseif W isa AbstractSciMLOperator
        update_coefficients!(W, uprev, p, t; gamma = h_fast)
    else
        jacobian2W!(W, f.mass_matrix, h_fast, cache.J)
    end
    integrator.stats.nw += 1
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
    ns = _extrapolation_sequence(alg.seq, order)

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
        # solves below reuses that factorization (`dolinsolve` without `A`).
        _mreil_update_W!(integrator, cache, h_fast)
        if !issuccess_W(W)
            OrdinaryDiffEqCore.set_EEst!(integrator, 2)
            return nothing
        end
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

                # Every column restarts from (uprev, t), so the first substep of the
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
                if !SciMLBase.successful_retcode(linres.retcode) &&
                        linres.retcode != SciMLBase.ReturnCode.Default
                    OrdinaryDiffEqCore.set_EEst!(integrator, 2)
                    return nothing
                end

                vecdz = _vec(linres.u)
                @.. broadcast = false vecTj = vecTj + vecdz
            end
        end
    end

    # Aitken–Neville Richardson extrapolation (in-place, reverse-row order)
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
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

    # Right endpoint derivative for the Hermite interpolant (k[2]).
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
    ns = _extrapolation_sequence(alg.seq, order)
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
            OrdinaryDiffEqCore.set_EEst!(integrator, 2)
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

                # Every column restarts from (uprev, t), so the first substep of the
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

    # Aitken–Neville Richardson extrapolation
    for k in 1:(order - 1)
        for j in order:-1:(k + 1)
            ratio = ns[j] / ns[j - k]
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

    # Right endpoint derivative for the Hermite interpolant (k[2]).
    integrator.fsallast = f.f1(integrator.u, p, t + dt) + f.f2(integrator.u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function initialize!(integrator, cache::MRABCache)
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

function initialize!(integrator, cache::MRABConstantCache)
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

# ── MRAB perform_step! (in-place, MutableCache) ───────────────────────────────
# Substep ℓ bootstraps with AB-min(ℓ, k) before the history reaches k samples.

function perform_step!(integrator, cache::MRABCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, k_slow, k_fast, F_history, tab) = cache
    β = tab.β
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = β[k]

    @.. broadcast = false u = uprev
    f.f2(k_slow, u, p, t)
    integrator.stats.nf2 += 1

    n_hist = 0
    for ℓ in 1:m
        t_fast = t + (ℓ - 1) * h
        f.f1(k_fast, u, p, t_fast)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        @.. broadcast = false F_history[min(n_hist + 1, k)] = k_fast + k_slow
        for j in min(n_hist + 1, k):-1:2
            F_history[j], F_history[j - 1] = F_history[j - 1], F_history[j]
        end
        n_hist = min(n_hist + 1, k)

        ks = min(ℓ, k)
        βs_use = β[ks]
        for i in 1:ks
            βi = βs_use[i]
            Fi = F_history[i]
            @.. broadcast = false u = u + h * βi * Fi
        end
    end

    # Error estimate (AB-k − AB-(k-1)) folds into one linear combination over F_history.
    return if integrator.opts.adaptive && k >= 2
        βs_low = β[k - 1]
        @.. broadcast = false tmp = h * (βs[1] - βs_low[1]) * F_history[1]
        for i in 2:(k - 1)
            @.. broadcast = false tmp = tmp + h * (βs[i] - βs_low[i]) * F_history[i]
        end
        @.. broadcast = false tmp = tmp + h * βs[k] * F_history[k]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

# ── MRAB perform_step! (out-of-place, ConstantCache) ──────────────────────────

@muladd function perform_step!(integrator, cache::MRABConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    β = cache.tab.β
    alg = unwrap_alg(integrator, false)
    k = alg.k
    m = alg.m
    h = dt / m
    βs = β[k]

    u_cur = uprev
    k_slow = f.f2(u_cur, p, t)
    integrator.stats.nf2 += 1

    F_history = Vector{typeof(k_slow)}(undef, k)
    n_hist = 0

    for ℓ in 1:m
        t_fast = t + (ℓ - 1) * h
        k_fast = f.f1(u_cur, p, t_fast)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        F = k_fast + k_slow
        n_hist = min(n_hist + 1, k)
        for j in n_hist:-1:2
            F_history[j] = F_history[j - 1]
        end
        F_history[1] = F

        ks = min(ℓ, k)
        βs_use = β[ks]
        Δu = zero(u_cur)
        for i in 1:ks
            Δu = @.. broadcast = false Δu + h * βs_use[i] * F_history[i]
        end
        u_cur = @.. broadcast = false u_cur + Δu
    end

    integrator.u = u_cur

    if integrator.opts.adaptive && k >= 2
        βs_low = β[k - 1]
        utilde = @.. broadcast = false h * (βs[1] - βs_low[1]) * F_history[1]
        for i in 2:(k - 1)
            utilde = @.. broadcast = false utilde + h * (βs[i] - βs_low[i]) * F_history[i]
        end
        utilde = @.. broadcast = false utilde + h * βs[k] * F_history[k]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

function initialize!(integrator, cache::MRIGARKCache)
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

function initialize!(integrator, cache::MRIGARKConstantCache)
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

# Fast-IVP rate at normalized τ: kbuf = Δcᵢ·dt·f1(w, t_phys) + dt·Σⱼ ωⱼ(τ)·fS[j],
# ωⱼ(τ) = w0[j] + w1[j]·τ (w1=nothing means constant coupling). IIP.
function _mrigark_rate!(kbuf, f1eval, f, w, p, t, dt, Δci, cprev, τ, w0, w1, fS, upto)
    t_phys = t + (cprev + τ * Δci) * dt
    f.f1(f1eval, w, p, t_phys)
    @.. broadcast = false kbuf = dt * Δci * f1eval
    for j in 1:upto
        ω = w1 === nothing ? w0[j] : w0[j] + w1[j] * τ
        iszero(ω) && continue
        @.. broadcast = false kbuf = kbuf + dt * ω * fS[j]
    end
    return kbuf
end

# Integrate one fast sub-stage from `vstart` over τ∈[0,1] into `out` (IIP).
# Δcᵢ=0 ⇒ pure slow quadrature (∫₀¹ω dτ = w0 + w1/2). q = inner explicit-RK order.
function _mrigark_substage!(
        out, vstart, cache, f, p, t, dt, Δci, cprev, m, q, w0, w1, fS, upto, stats
    )
    (; v, vtmp, f1eval, kk) = cache
    if iszero(Δci)
        @.. broadcast = false out = vstart
        for j in 1:upto
            ω = w1 === nothing ? w0[j] : w0[j] + w1[j] / 2
            iszero(ω) && continue
            @.. broadcast = false out = out + dt * ω * fS[j]
        end
        return out
    end
    h = 1 / m
    @.. broadcast = false v = vstart
    for ks in 1:m
        τ0 = (ks - 1) * h
        _mrigark_rate!(kk[1], f1eval, f, v, p, t, dt, Δci, cprev, τ0, w0, w1, fS, upto)
        @.. broadcast = false vtmp = v + (h / 2) * kk[1]
        _mrigark_rate!(kk[2], f1eval, f, vtmp, p, t, dt, Δci, cprev, τ0 + h / 2, w0, w1, fS, upto)
        if q == 2
            @.. broadcast = false v = v + h * kk[2]
        elseif q == 3
            @.. broadcast = false vtmp = v - h * kk[1] + 2 * h * kk[2]
            _mrigark_rate!(kk[3], f1eval, f, vtmp, p, t, dt, Δci, cprev, τ0 + h, w0, w1, fS, upto)
            @.. broadcast = false v = v + (h / 6) * (kk[1] + 4 * kk[2] + kk[3])
        else
            @.. broadcast = false vtmp = v + (h / 2) * kk[2]
            _mrigark_rate!(kk[3], f1eval, f, vtmp, p, t, dt, Δci, cprev, τ0 + h / 2, w0, w1, fS, upto)
            @.. broadcast = false vtmp = v + h * kk[3]
            _mrigark_rate!(kk[4], f1eval, f, vtmp, p, t, dt, Δci, cprev, τ0 + h, w0, w1, fS, upto)
            @.. broadcast = false v = v + (h / 6) * (kk[1] + 2 * kk[2] + 2 * kk[3] + kk[4])
        end
        OrdinaryDiffEqCore.increment_nf!(stats, q)
    end
    @.. broadcast = false out = v
    return out
end

function perform_step!(integrator, cache::MRIGARKCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, z, fS, zemb, tab) = cache
    (; Δc, W0, W1, Wemb0, Wemb1, q) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Δc)
    stats = integrator.stats

    @.. broadcast = false z[1] = uprev
    cprev = zero(eltype(Δc))
    for i in 1:s
        f.f2(fS[i], z[i], p, t + cprev * dt)
        stats.nf2 += 1
        _mrigark_substage!(
            z[i + 1], z[i], cache, f, p, t, dt, Δc[i], cprev, m, q,
            view(W0, i, :), view(W1, i, :), fS, i, stats
        )
        cprev += Δc[i]
    end
    @.. broadcast = false u = z[s + 1]

    return if integrator.opts.adaptive
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
end

# Fast-IVP rate, out-of-place.
function _mrigark_rate(f, w, p, t, dt, Δci, cprev, τ, w0, w1, fS, upto)
    t_phys = t + (cprev + τ * Δci) * dt
    ffast = f.f1(w, p, t_phys)
    g = @.. broadcast = false dt * Δci * ffast
    for j in 1:upto
        ω = w1 === nothing ? w0[j] : w0[j] + w1[j] * τ
        iszero(ω) && continue
        g = @.. broadcast = false g + dt * ω * fS[j]
    end
    return g
end

function _mrigark_substage(vstart, f, p, t, dt, Δci, cprev, m, q, w0, w1, fS, upto, stats)
    if iszero(Δci)
        v = vstart
        for j in 1:upto
            ω = w1 === nothing ? w0[j] : w0[j] + w1[j] / 2
            iszero(ω) && continue
            v = @.. broadcast = false v + dt * ω * fS[j]
        end
        return v
    end
    h = 1 / m
    v = vstart
    for ks in 1:m
        τ0 = (ks - 1) * h
        k1 = _mrigark_rate(f, v, p, t, dt, Δci, cprev, τ0, w0, w1, fS, upto)
        vmid = @.. broadcast = false v + (h / 2) * k1
        k2 = _mrigark_rate(f, vmid, p, t, dt, Δci, cprev, τ0 + h / 2, w0, w1, fS, upto)
        if q == 2
            v = @.. broadcast = false v + h * k2
        elseif q == 3
            vmid2 = @.. broadcast = false v - h * k1 + 2 * h * k2
            k3 = _mrigark_rate(f, vmid2, p, t, dt, Δci, cprev, τ0 + h, w0, w1, fS, upto)
            v = @.. broadcast = false v + (h / 6) * (k1 + 4 * k2 + k3)
        else
            vmid2 = @.. broadcast = false v + (h / 2) * k2
            k3 = _mrigark_rate(f, vmid2, p, t, dt, Δci, cprev, τ0 + h / 2, w0, w1, fS, upto)
            vmid3 = @.. broadcast = false v + h * k3
            k4 = _mrigark_rate(f, vmid3, p, t, dt, Δci, cprev, τ0 + h, w0, w1, fS, upto)
            v = @.. broadcast = false v + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        end
        OrdinaryDiffEqCore.increment_nf!(stats, q)
    end
    return v
end

@muladd function perform_step!(
        integrator, cache::MRIGARKConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; Δc, W0, W1, Wemb0, Wemb1, q) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Δc)
    stats = integrator.stats

    z = Vector{typeof(uprev)}(undef, s + 1)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)
    z[1] = uprev
    cprev = zero(eltype(Δc))
    for i in 1:s
        fS[i] = f.f2(z[i], p, t + cprev * dt)
        stats.nf2 += 1
        z[i + 1] = _mrigark_substage(
            z[i], f, p, t, dt, Δc[i], cprev, m, q,
            view(W0, i, :), view(W1, i, :), fS, i, stats
        )
        cprev += Δc[i]
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

    return if integrator.opts.adaptive
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
end

@muladd function perform_step!(
        integrator, cache::MRIGARKImplicitConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; nlsolver, tab) = cache
    (; Δc, W0, W1, Wemb0, Wemb1, γ0, q) = tab
    alg = unwrap_alg(integrator, true)
    m = alg.m
    s = length(Δc)
    stats = integrator.stats

    markfirststage!(nlsolver)

    z = Vector{typeof(uprev)}(undef, s + 1)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)
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
end

# Stage 1 is identity (Y_1 = u_n); inner ODE active only for i ≥ 2.

function initialize!(integrator, cache::MISCache)
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

function initialize!(integrator, cache::MISConstantCache)
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

function perform_step!(integrator, cache::MISCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, atmp, v, offset, k_fast, Y, fS, tab) = cache
    (; α, β, γ, d, c, ctilde) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(d)

    @.. broadcast = false Y[1] = uprev
    f.f2(fS[1], uprev, p, t)
    integrator.stats.nf2 += 1

    for i in 2:s
        d_i = d[i]
        inv_di = 1 / d_i

        @.. broadcast = false v = uprev
        for j in 1:(i - 1)
            if !iszero(α[i, j])
                @.. broadcast = false v = v + α[i, j] * (Y[j] - uprev)
            end
        end

        @.. broadcast = false offset = zero(offset)
        for j in 1:(i - 1)
            if !iszero(γ[i, j])
                @.. broadcast = false offset = offset + (γ[i, j] * inv_di / dt) * (Y[j] - uprev)
            end
            if !iszero(β[i, j])
                @.. broadcast = false offset = offset + (β[i, j] * inv_di) * fS[j]
            end
        end

        M_inner = max(1, ceil(Int, m * d_i))
        h_inner = d_i * dt / M_inner
        slope = (c[i] - ctilde[i]) * inv_di
        for k_step in 1:M_inner
            τ = (k_step - 1) * h_inner
            t_a = t + ctilde[i] * dt + slope * τ
            t_b = t + ctilde[i] * dt + slope * (τ + h_inner / 2)
            f.f1(k_fast, v, p, t_a)
            @.. broadcast = false tmp = v + (h_inner / 2) * (offset + k_fast)
            f.f1(k_fast, tmp, p, t_b)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            @.. broadcast = false v = v + h_inner * (offset + k_fast)
        end

        @.. broadcast = false Y[i] = v
        f.f2(fS[i], Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1
    end

    @.. broadcast = false u = Y[s]

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = Y[s] - Y[s - 1]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(integrator, cache::MISConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    (; α, β, γ, d, c, ctilde) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(d)

    Y = Vector{typeof(uprev)}(undef, s)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)

    Y[1] = uprev
    fS[1] = f.f2(uprev, p, t)
    integrator.stats.nf2 += 1

    for i in 2:s
        d_i = d[i]
        inv_di = 1 / d_i

        v_state = uprev
        for j in 1:(i - 1)
            if !iszero(α[i, j])
                v_state = @.. broadcast = false v_state + α[i, j] * (Y[j] - uprev)
            end
        end

        offset = zero(fS[1])
        for j in 1:(i - 1)
            if !iszero(γ[i, j])
                offset = @.. broadcast = false offset + (γ[i, j] * inv_di / dt) * (Y[j] - uprev)
            end
            if !iszero(β[i, j])
                offset = @.. broadcast = false offset + (β[i, j] * inv_di) * fS[j]
            end
        end

        M_inner = max(1, ceil(Int, m * d_i))
        h_inner = d_i * dt / M_inner
        slope = (c[i] - ctilde[i]) * inv_di
        for k_step in 1:M_inner
            τ = (k_step - 1) * h_inner
            t_a = t + ctilde[i] * dt + slope * τ
            t_b = t + ctilde[i] * dt + slope * (τ + h_inner / 2)
            k_fast = f.f1(v_state, p, t_a)
            v_mid = @.. broadcast = false v_state + (h_inner / 2) * (offset + k_fast)
            k_fast = f.f1(v_mid, p, t_b)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            v_state = @.. broadcast = false v_state + h_inner * (offset + k_fast)
        end

        Y[i] = v_state
        fS[i] = f.f2(Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1
    end

    integrator.u = Y[s]

    if integrator.opts.adaptive
        utilde = @.. broadcast = false Y[s] - Y[s - 1]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

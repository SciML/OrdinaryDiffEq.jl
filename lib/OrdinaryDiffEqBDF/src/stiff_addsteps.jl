####################################################################
# BDF _ode_addsteps!: rebuild interpolation data after callbacks
#
# _ode_addsteps! is called in two contexts:
#   1. After callbacks (always_calc_begin = true): k needs rebuilding
#      from the cache's current history for the truncated step.
#   2. During post-solve interpolation (always_calc_begin = false):
#      k was already saved correctly by perform_step!, and the cache
#      state corresponds to the final step, not the step being
#      interpolated. In this case we must be a no-op.
####################################################################

####################################################################
# QNDF: Rebuild backward differences for interpolation
#
# Layout: k[1..max_order] = backward differences D[j] on the grid of the
# completed step, of length h = cache.dtprev:
#     p(Θ) = y_h + Σ_{m=1}^{order} φ_m(Θ-1) D[m],   Θ = (t - tprev) / h,
# with y_h = uprev + D[1] the end of the completed step and
# φ_m(σ) = σ (σ+1) ⋯ (σ+m-1) / m!.
#
# A callback truncates the step to dt = θ h (θ = dt / h). The interpolant
# over the truncated step is the same polynomial in Θ' = (t - tprev) / dt,
# p(θ Θ'), whose backward differences on the grid of spacing dt ending at
# the new endpoint are
#     k[j] = Σ_{i=0}^{j} (-1)^i C(j, i) p(θ (1 - i))
#          = Σ_{m=j}^{order} T(θ, m, j) D[m],
#     T(θ, m, j) = Σ_{i=0}^{j} (-1)^i C(j, i) φ_m(θ (1 - i) - 1)
# (the j-th difference of a polynomial of degree m < j vanishes; T = I at
# θ = 1). The first difference is k[1] = u - uprev, which is p(θ) - p(0)
# unless the callback modified u: the interpolant then passes through the
# modified endpoint.
####################################################################

# φ_m(σ) = σ (σ+1) ⋯ (σ+m-1) / m!
function _qndf_phi(σ, m)
    φ = one(σ)
    for j in 1:m
        φ *= (σ + j - 1) / j
    end
    return φ
end

# The coefficient of D[m] in the j-th backward difference of the step truncated to
# θ of the completed one, T(θ, m, j)
function _qndf_rebase_coefficient(θ, m, j)
    c = zero(θ)
    for i in 0:j
        c += (-1)^i * binomial(j, i) * _qndf_phi(θ * (1 - i) - 1, m)
    end
    return c
end

# QNDF ConstantCache: out-of-place k entries
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::QNDFConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    always_calc_begin || return nothing
    (; D, order, dtprev) = cache
    θ = dt / dtprev
    for j in eachindex(k)
        if j == 1
            # First difference must match new endpoints: k[1] = u - uprev
            k[j] = u isa Number ? (u - uprev) : @.. u - uprev
        elseif j <= order
            # The completed step's polynomial on the truncated step
            kj = zero(u)
            for m in j:order
                c = _qndf_rebase_coefficient(θ, m, j)
                kj = u isa Number ? kj + c * D[m] : @.. kj + c * D[m]
            end
            k[j] = kj
        else
            k[j] = zero(u)
        end
    end
    return nothing
end

# QNDF Cache: in-place k entries (pre-allocated arrays)
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::QNDFCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    always_calc_begin || return nothing
    (; D, order, dtprev) = cache
    θ = dt / dtprev
    @.. broadcast = false k[1] = u - uprev
    for j in 2:length(k)
        fill!(k[j], zero(eltype(u)))
        j <= order || continue
        # The completed step's polynomial on the truncated step
        for m in j:order
            c = _qndf_rebase_coefficient(θ, m, j)
            @.. broadcast = false k[j] += c * D[m]
        end
    end
    return nothing
end

####################################################################
# FBDF / DFBDF: Rebuild Chebyshev-resampled interpolation data
#
# Layout: k has max_order+1 entries.
#   k[1..n] = solution values at fixed Chebyshev reference nodes
#   k[n+1..end] = zero
#
# After a callback truncates the step to dt (= t_event - t_prev)
# and modifies u, we rebuild k by:
#   1. Assembling the original Lagrange interpolation data:
#      values[1] = u at Θ=1, values[1+j] = u_history[j]
#      thetas[1] = 1, thetas[1+j] = (ts[j] - t) / dt
#   2. Resampling at fixed Chebyshev nodes via _resample_at_chebyshev!
####################################################################

# FBDF/DFBDF ConstantCache: out-of-place k entries
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{FBDFConstantCache, DFBDFConstantCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    always_calc_begin || return nothing
    (; u_history, ts, order) = cache
    n = order + 1

    thetas = Vector{typeof(t)}(undef, n)
    thetas[1] = one(t)
    for j in 1:order
        thetas[1 + j] = (ts[j] - t) / dt
    end

    values = Vector{typeof(u)}(undef, n)
    values[1] = u isa Number ? u : copy(u)
    for j in 1:order
        values[1 + j] = u isa Number ? u_history[j] : copy(u_history[j])
    end

    _resample_at_chebyshev!(k, values, thetas, n)
    for j in (n + 1):length(k)
        k[j] = zero(u)
    end
    return nothing
end

# FBDF/DFBDF Cache: in-place k entries (pre-allocated arrays)
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{FBDFCache, DFBDFCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    always_calc_begin || return nothing
    (; u_history, ts, order, equi_ts) = cache
    max_order = length(k) - 1
    n = order + 1

    # Compute thetas
    equi_ts[1] = one(eltype(equi_ts))
    for j in 1:order
        equi_ts[1 + j] = (ts[j] - t) / dt
    end

    # Resample at Chebyshev nodes directly from u and u_history
    _resample_at_chebyshev_direct_iip!(k, u, u_history, equi_ts, n)
    for j in (n + 1):(max_order + 1)
        fill!(k[j], zero(eltype(u)))
    end
    return nothing
end

####################################################################
# NordsieckBDF / DNordsieckBDF: rebase the Nordsieck polynomial
#
# k[j] = zn[j] are the columns about the endpoint of the step h = hscale:
#     p(t) = Σ_j k[j] s^(j-1),   s = (t - (tprev + h)) / h.
# A callback shortening the step to dt evaluates the interpolant with
# s' = (t - (tprev + dt)) / dt, and s = a + b s' with b = dt/h, a = b - 1.
# Rebuild from zn (not k) so repeated calls stay idempotent.
#
# While `cache.predicted` holds (post-perform_step!, pre-acceptance) zn is the
# uncommitted trial state, and a forced rebuild yields the corrected trial
# polynomial zn[j] + l[j]*acor — what nordsieck_complete! commits on
# acceptance. DelayDiffEq's advance/update of the history integrator relies on
# this; evaluation-only callers (always_calc_begin = force_calc_end = false)
# keep whatever k already holds.
####################################################################

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{NordsieckBDFConstantCache, DNordsieckBDFConstantCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    (; zn, order, hscale) = cache
    iszero(hscale) && return nothing
    trial = cache.predicted
    (always_calc_begin || (trial && force_calc_end)) || return nothing
    b = dt / hscale
    a = b - one(b)
    n = length(k)
    # constant caches alias k[j] to zn[j]; never mutate those in place
    for j in 1:n
        k[j] = if j - 1 > order
            zero(u)
        elseif trial
            @.. zn[j] + cache.l[j] * cache.acor
        else
            zn[j]
        end
    end
    for i in 1:(n - 1)
        for j in (n - 1):-1:i
            k[j] = @.. k[j] + a * k[j + 1]
        end
    end
    scale = b
    for j in 2:n
        k[j] = @.. scale * k[j]
        scale *= b
    end
    return nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{NordsieckBDFCache, DNordsieckBDFCache},
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    (; zn, order, hscale) = cache
    iszero(hscale) && return nothing
    trial = cache.predicted
    (always_calc_begin || (trial && force_calc_end)) || return nothing
    b = dt / hscale
    a = b - one(b)
    n = length(k)
    for j in 1:n
        if j - 1 <= order
            if trial
                kj = k[j]
                @.. broadcast = false kj = zn[j] + cache.l[j] * cache.acor
            else
                copyto!(k[j], zn[j])
            end
        else
            fill!(k[j], zero(eltype(u)))
        end
    end
    for i in 1:(n - 1)
        for j in (n - 1):-1:i
            kj, kj1 = k[j], k[j + 1]
            @.. broadcast = false kj = kj + a * kj1
        end
    end
    scale = b
    for j in 2:n
        kj = k[j]
        @.. broadcast = false kj = scale * kj
        scale *= b
    end
    return nothing
end

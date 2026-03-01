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
# Layout: k[1..max_order] = backward differences D[j]
# Interpolation: p(Θ) = y₁ + Σ φ_j(Θ-1) * k[j]
# where y₁ = u (step endpoint).
#
# After a callback truncates the step and modifies u, we rebuild k
# using the cache's D matrix for higher-order accuracy. The first
# difference k[1] is updated to u - uprev so that the interpolant
# passes through the correct endpoints; higher differences D[j≥2]
# are preserved from the completed step.
####################################################################

# QNDF ConstantCache: out-of-place k entries
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::QNDFConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    always_calc_begin || return nothing
    (; D, order) = cache
    for j in eachindex(k)
        if j == 1
            # First difference must match new endpoints: k[1] = u - uprev
            k[j] = u isa Number ? (u - uprev) : @.. u - uprev
        elseif j <= order
            # Preserve higher-order differences from the completed step
            k[j] = D[j]
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
    (; D, order) = cache
    @.. broadcast = false k[1] = u - uprev
    for j in 2:length(k)
        if j <= order
            copyto!(k[j], D[j])
        else
            fill!(k[j], zero(eltype(u)))
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

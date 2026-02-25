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
# FBDF / DFBDF: Rebuild Lagrange interpolation nodes
#
# Layout: k has 2*half entries where half = length(k)÷2.
#   k[1..half]       = solution values at nodes
#   k[half+1..2*half] = corresponding Θ positions
#
# After a callback truncates the step to dt (= t_event - t_prev)
# and modifies u, we rebuild k using the cache's u_history and ts:
#   k[1] = u (post-callback) at Θ = 1
#   k[1+j] = u_history[j] (past solutions, still valid)
#   k[half+1] = 1
#   k[half+1+j] = (ts[j] - t) / dt  (recomputed for truncated dt)
#
# This preserves high-order Lagrange interpolation accuracy through
# the actual historical solution values with correct Θ normalization.
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
    half = length(k) ÷ 2
    if u isa Number
        k[1] = u
        k[half + 1] = one(t)
        for j in 1:(half - 1)
            if j <= order
                k[1 + j] = u_history[j]
                k[half + 1 + j] = (ts[j] - t) / dt
            else
                k[1 + j] = zero(u)
                k[half + 1 + j] = zero(t)
            end
        end
    else
        k[1] = copy(u)
        k[half + 1] = fill(one(t), size(u))
        for j in 1:(half - 1)
            if j <= order
                k[1 + j] = copy(u_history[j])
                k[half + 1 + j] = fill((ts[j] - t) / dt, size(u))
            else
                k[1 + j] = zero(u)
                k[half + 1 + j] = fill(zero(t), size(u))
            end
        end
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
    (; u_history, ts, order) = cache
    half = length(k) ÷ 2
    @.. broadcast = false k[1] = u
    fill!(k[half + 1], one(eltype(u)))
    for j in 1:(half - 1)
        if j <= order
            copyto!(k[1 + j], u_history[j])
            fill!(k[half + 1 + j], (ts[j] - t) / dt)
        else
            fill!(k[1 + j], zero(eltype(u)))
            fill!(k[half + 1 + j], zero(eltype(u)))
        end
    end
    return nothing
end

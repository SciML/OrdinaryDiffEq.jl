@concrete mutable struct NewmarkDiscretizationCache
    # Eval
    f
    t
    p
    # Discretization params (αm = αf = 0 recovers Newmark)
    dt
    β
    γ
    αm
    αf
    aₙ
    vₙ
    uₙ
    # # Buffers be compatible with the chosen AD mode
    atmp
    vₙ₊₁
    uₙ₊₁
end

# This is derived from the idea stated in Nonlinear Finite Elements by Peter Wriggers, Ch 6.1.2 .
#
# Let us introduce the notation v = u' and a = u'' = v' such that we write the ODE problem as Ma = f(v,u,t).
# For the time discretization we assume that:
#   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 aₙ₊ₐ₁
#   vₙ₊₁ = vₙ + Δtₙ aₙ₊ₐ₂
# with a₁ = 1-2β and a₂ = 1-γ, such that
#   uₙ₊₁ = uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁]
#   vₙ₊₁ = vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁]
#
# For Generalized-α the equations of motion are evaluated at interpolated states:
#   M·aₙ₊αm = f(uₙ₊αf, vₙ₊αf, tₙ₊αf)
# Setting αm = αf = 0 recovers Newmark exactly
#
# For the Newton method the effective Jacobian is:
#   (1-αm)·M - (1-αf)·(Δtₙ²β ∂fᵤ + Δtₙγ ∂fᵥ) = 0

# in place variant
@muladd function discretized_residual!(
        residual, aₙ₊₁, cache::NewmarkDiscretizationCache
    )
    (; f, dt, t, p) = cache
    (; β, γ, αm, αf, aₙ, vₙ, uₙ) = cache

    atmp = get_tmp(cache.atmp, aₙ₊₁)
    vₙ₊₁ = get_tmp(cache.vₙ₊₁, aₙ₊₁)
    uₙ₊₁ = get_tmp(cache.uₙ₊₁, aₙ₊₁)

    # standard newmark update for full step quantities
    @.. uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    # interpolates to αf for the state evaluation blending
    @.. uₙ₊₁ = (1 - αf) * uₙ₊₁ + αf * uₙ
    @.. vₙ₊₁ = (1 - αf) * vₙ₊₁ + αf * vₙ
    tₙ₊αf = t + (1 - αf) * dt

    f.f1(atmp, vₙ₊₁, uₙ₊₁, p, tₙ₊αf)
    M = f.mass_matrix

    mul!(residual, M, (1 - αm) * aₙ₊₁ + αm * aₙ)
    @.. residual = residual - atmp

    return nothing
end

# out of-place variant
@muladd function discretized_residual(aₙ₊₁, cache::NewmarkDiscretizationCache)
    (; f, dt, t, p) = cache
    (; β, γ, αm, αf, aₙ, vₙ, uₙ) = cache

    uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    uₙ₊αf = (1 - αf) * uₙ₊₁ + αf * uₙ
    vₙ₊αf = (1 - αf) * vₙ₊₁ + αf * vₙ
    tₙ₊αf = t + (1 - αf) * dt

    aₙ₊αm = (1 - αm) * aₙ₊₁ + αm * aₙ

    M = f.mass_matrix
    return M * aₙ₊αm - f.f1(vₙ₊αf, uₙ₊αf, p, tₙ₊αf)
end

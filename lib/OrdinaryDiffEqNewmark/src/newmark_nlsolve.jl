@concrete mutable struct NewmarkDiscretizationCache
    # Eval
    f
    t
    p
    # Newmark discretization params
    dt
    β
    γ
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
# This allows us to reduce the implicit discretization to have only aₙ₊₁ as the unknown:
#   Maₙ₊₁ = f(vₙ₊₁(aₙ₊₁), uₙ₊₁(aₙ₊₁), tₙ₊₁)
#         = f(vₙ + Δtₙ [(1-γ)aₙ + γaₙ₊₁], uₙ + Δtₙ vₙ + Δtₙ²/2 [(1-2β)aₙ + 2βaₙ₊₁], tₙ₊₁)
# Such that we have to solve the nonlinear problem
#   Maₙ₊₁ - f(vₙ₊₁(aₙ₊₁), uₙ₊₁(aₙ₊₁), tₙ₊₁)  = 0
# for aₙ₊₁'' in each time step.

# For the Newton method the linearization becomes
#   M - (dₐuₙ₊₁ ∂fᵤ + dₐvₙ₊₁ ∂fᵥ) = 0
#   M - (Δtₙ²β  ∂fᵤ +  Δtₙγ  ∂fᵥ) = 0

# Inplace variant
@muladd function newmark_discretized_residual!(
        residual, aₙ₊₁, p_newmark::NewmarkDiscretizationCache
    )
    (; f, dt, t, p) = p_newmark
    (; γ, β, aₙ, vₙ, uₙ) = p_newmark

    atmp = get_tmp(p_newmark.atmp, aₙ₊₁)
    vₙ₊₁ = get_tmp(p_newmark.vₙ₊₁, aₙ₊₁)
    uₙ₊₁ = get_tmp(p_newmark.uₙ₊₁, aₙ₊₁)

    @.. uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    @.. vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    # This temporary variable is also not compatible with AD.
    f.f1(atmp, vₙ₊₁, uₙ₊₁, p, t)
    M = f.mass_matrix

    mul!(residual, M, aₙ₊₁)
    @.. residual = residual - atmp

    return nothing
end

# Out of place variant
@muladd function newmark_discretized_residual(aₙ₊₁, p_newmark::NewmarkDiscretizationCache)
    (; f, dt, t, p) = p_newmark
    (; γ, β, aₙ, vₙ, uₙ) = p_newmark

    uₙ₊₁ = uₙ + dt * vₙ + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    vₙ₊₁ = vₙ + dt * ((1 - γ) * aₙ + γ * aₙ₊₁)

    M = f.mass_matrix
    return M * aₙ₊₁ - f.f1(vₙ₊₁, uₙ₊₁, p, t)
end

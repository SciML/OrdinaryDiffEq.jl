# ============================================================================
# Runtime Tableauization for Generic Rosenbrock Methods
# ============================================================================
#
# This file implements a runtime-based approach for Rosenbrock methods,
# similar to how ExplicitRK handles tableauization. Instead of using
# compile-time code generation (macros), this approach:
#
# 1. Stores tableau coefficients (A, C, gamma, c, d, H) in arrays
# 2. Uses runtime loops to compute stages
# 3. Supports any tableau implementing the RodasTableau interface
#
# Benefits:
# - Faster compile times (no macro expansion)
# - Easier to add new methods (just define a tableau)
# - Simpler code maintenance
#
# Trade-offs:
# - Slightly slower runtime due to loop overhead (usually negligible)
# ============================================================================

# ============================================================================
# Runtime Cache Structures
# ============================================================================

"""
    GenericRosenbrockRuntimeConstantCache

Constant cache for runtime GenericRosenbrock solver (out-of-place problems).
Stores the tableau and pre-allocated arrays for stage values.
"""
struct GenericRosenbrockRuntimeConstantCache{
    TabType, TF, UF, JType, WType, F, AD
} <: RosenbrockConstantCache
    tab::TabType       # RodasTableau with A, C, gamma, c, d, H
    tf::TF             # Time derivative wrapper
    uf::UF             # State derivative wrapper
    J::JType           # Jacobian matrix
    W::WType           # W = I - dt*gamma*J
    linsolve::F        # Linear solver
    autodiff::AD       # Autodiff choice
    interp_order::Int  # Number of interpolation polynomial coefficients
end

# Mutable cache for runtime GenericRosenbrock solver (in-place problems).
# Stores pre-allocated working arrays for stages and intermediate values.
@cache mutable struct GenericRosenbrockRuntimeCache{
    uType, rateType, TabType, uNoUnitsType, JType, WType,
    TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter
} <: GenericRosenbrockMutableCache
    u::uType
    uprev::uType
    dense::Vector{rateType}           # Dense output coefficients (H-transformed)
    ks::Vector{rateType}              # Stage values k1, k2, ..., kn
    du::rateType
    du1::rateType
    du2::rateType
    f_tmp::rateType               # Temporary for function evaluations
    fsalfirst::rateType
    fsallast::rateType
    dT::rateType                  # Time derivative
    J::JType
    W::WType
    tmp::rateType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    tab::TabType
    tf::TFType
    uf::UFType
    linsolve_tmp::rateType
    linsolve::F
    jac_config::JCType
    grad_config::GCType
    reltol::RTolType
    alg::A
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    interp_order::Int
end

# ============================================================================
# Cache Construction
# ============================================================================

"""
Build the runtime cache for in-place (mutating) problems.
"""
function alg_cache(
    alg::GenericRosenbrock, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}, verbose
) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = alg.tableau
    stages = size(tab.A, 1)
    interp_order = size(tab.H, 1)

    # Pre-allocate dense output arrays
    dense = [zero(rate_prototype) for _ in 1:interp_order]

    # Pre-allocate stage vectors
    ks = [zero(rate_prototype) for _ in 1:stages]

    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    f_tmp = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    linsolve_tmp = zero(rate_prototype)

    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(true))

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve = nothing
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)
    grad_config = build_grad_config(alg, f, tf, du1, t)

    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))
    Pl, Pr = wrapprecs(
        alg.precs(W, nothing, u, p, t, nothing, nothing, nothing, nothing)..., weight, tmp
    )
    linsolve = init(linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity)

    GenericRosenbrockRuntimeCache(
        u, uprev, dense, ks, du, du1, du2, f_tmp, fsalfirst, fsallast, dT,
        J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp, linsolve,
        jac_config, grad_config, reltol, alg, alg.step_limiter!, alg.stage_limiter!,
        interp_order
    )
end

"""
Build the runtime cache for out-of-place (non-mutating) problems.
"""
function alg_cache(
    alg::GenericRosenbrock, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}, verbose
) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = alg.tableau

    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)

    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))

    GenericRosenbrockRuntimeConstantCache(
        tab, tf, uf, J, W, nothing, alg.autodiff, size(tab.H, 1)
    )
end

# ============================================================================
# Initialize
# ============================================================================

function initialize!(integrator, cache::GenericRosenbrockRuntimeCache)
    integrator.kshortsize = cache.interp_order
    resize!(integrator.k, integrator.kshortsize)
    for i in 1:integrator.kshortsize
        integrator.k[i] = cache.dense[i]
    end
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::GenericRosenbrockRuntimeConstantCache)
    integrator.kshortsize = cache.interp_order
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    for i in 1:integrator.kshortsize
        integrator.k[i] = zero(integrator.u)
    end
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
end

# ============================================================================
# Runtime perform_step! (In-Place Version)
# ============================================================================

"""
    perform_step!(integrator, cache::GenericRosenbrockRuntimeCache, repeat_step=false)

Generic runtime implementation of Rosenbrock method for in-place problems.
Uses loops over the tableau instead of compile-time generated code.

The Rosenbrock method computes:
    W = I - dt*γ*J
    For i = 1, ..., s:
        W * k_i = f(u + Σⱼ a_ij*k_j) + dt*d_i*(∂f/∂t) + Σⱼ (C_ij/dt)*k_j
    u_new = u + Σᵢ b_i*k_i  (where b_i is the last row of A for stiffly accurate methods)
"""
@muladd function perform_step!(
    integrator, cache::GenericRosenbrockRuntimeCache, repeat_step = false
)
    (; t, dt, uprev, u, f, p) = integrator
    (; ks, du, du1, du2, f_tmp, fsalfirst, fsallast, dT, J, W, tmp, atmp, weight) = cache
    (; tab, tf, uf, linsolve_tmp, linsolve, step_limiter!, stage_limiter!) = cache

    A = tab.A
    C = tab.C
    gamma = tab.gamma
    c = tab.c
    d = tab.d

    stages = size(A, 1)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtgamma = dt * gamma

    # Precompute dt*C_ij and dt*d_i
    dtC = C ./ dt
    dtd = dt .* d

    f(fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute Jacobian and W matrix
    calc_rosenbrock_differentiation!(integrator, cache, dtgamma, dtgamma, repeat_step)

    # Compute weight for linear solver
    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    # ---- Stage 1 ----
    @.. broadcast=false linsolve_tmp = fsalfirst + dtd[1] * dT

    if repeat_step
        linres = dolinsolve(
            integrator, linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    else
        linres = dolinsolve(
            integrator, linsolve; A = W, b = _vec(linsolve_tmp),
            du = fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    end
    @.. broadcast=false $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    # ---- Stages 2 to s ----
    for i in 2:stages
        # Compute u = uprev + Σⱼ a_ij*k_j
        @.. broadcast=false u = uprev
        for j in 1:(i - 1)
            @.. broadcast=false u = u + A[i, j] * ks[j]
        end

        stage_limiter!(u, integrator, p, t + c[i] * dt)

        # Compute f(u, t + c_i*dt)
        f(du, u, p, t + c[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        # Build C terms: Σⱼ (C_ij/dt)*k_j
        @.. broadcast=false du1 = zero(eltype(du1))
        if mass_matrix === I
            for j in 1:(i - 1)
                @.. broadcast=false du1 = du1 + dtC[i, j] * ks[j]
            end
        else
            for j in 1:(i - 1)
                @.. broadcast=false du1 = du1 + dtC[i, j] * ks[j]
            end
            mul!(_vec(du2), mass_matrix, _vec(du1))
            @.. broadcast=false du1 = du2
        end
        @.. broadcast=false linsolve_tmp = du + dtd[i] * dT + du1

        # Solve W * k[i] = -linsolve_tmp
        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        @.. broadcast=false $(_vec(ks[i])) = -linres.u
        integrator.stats.nsolve += 1
    end

    # ---- Final solution ----
    # For methods with dense output stages (e.g., Rodas6P with 19 stages but
    # solution at stage 16), we must use the correct solution stage, not the last.
    sol_stage = cache.alg.num_solution_stages

    if sol_stage < stages
        # Recompute u from uprev using the solution stage row of A
        @.. broadcast=false du = ks[sol_stage]
        @.. broadcast=false u = uprev
        for j in 1:(sol_stage - 1)
            @.. broadcast=false u = u + A[sol_stage, j] * ks[j]
        end
        @.. broadcast=false u = u + ks[sol_stage]
    else
        # u already contains uprev + Σ A[s,j]*k[j] for j=1:s-1 from the last stage loop
        @.. broadcast=false du = ks[stages]
        @.. broadcast=false u = u + ks[stages]
    end

    step_limiter!(u, integrator, p, t + dt)

    # Compute fsallast for FSAL
    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # ---- Error estimation ----
    if integrator.opts.adaptive
        # For Rodas methods, the error estimate is the solution stage value
        calculate_residuals!(
            atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    # ---- Dense output coefficients ----
    if integrator.opts.calck
        H = tab.H
        for j in eachindex(integrator.k)
            integrator.k[j] .= 0
        end
        for i in eachindex(ks)
            for j in eachindex(integrator.k)
                @.. integrator.k[j] += H[j, i] * ks[i]
            end
        end
    end

    cache.linsolve = linres.cache
end

# ============================================================================
# Runtime perform_step! (Out-of-Place / Constant Cache Version)
# ============================================================================

@muladd function perform_step!(
    integrator, cache::GenericRosenbrockRuntimeConstantCache, repeat_step = false
)
    (; t, dt, uprev, u, f, p) = integrator
    (; tab, tf, uf) = cache

    A = tab.A
    C = tab.C
    gamma = tab.gamma
    c = tab.c
    d = tab.d

    stages = size(A, 1)
    mass_matrix = integrator.f.mass_matrix

    dtgamma = dt * gamma

    # Precompute dt-scaled coefficients
    dtC = C ./ dt
    dtd = dt .* d

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    # Pre-allocate ks array
    ks = Vector{typeof(integrator.fsalfirst)}(undef, stages)

    # ---- Stage 1 ----
    # Use f(uprev, p, t) directly like the reference implementation
    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    linsolve_tmp = @.. du + dtd[1] * dT
    ks[1] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    # ---- Stages 2 to s ----
    # u will be preserved from the last iteration for final solution
    u = uprev
    for i in 2:stages
        # Compute u = uprev + Σⱼ a_ij*k_j
        u = uprev
        for j in 1:(i - 1)
            u = @.. u + A[i, j] * ks[j]
        end

        # Compute f(u, t + c_i*dt)
        du = f(u, p, t + c[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        # Build RHS: f + dt*d_i*dT + Σⱼ (C_ij/dt)*k_j
        linsolve_tmp = zero(du)
        if mass_matrix === I
            for j in 1:(i - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[i, j] * ks[j]
            end
        else
            for j in 1:(i - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[i, j] * ks[j]
            end
            linsolve_tmp = mass_matrix * linsolve_tmp
        end
        linsolve_tmp = @.. du + dtd[i] * dT + linsolve_tmp

        ks[i] = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
        integrator.stats.nsolve += 1
    end

    # ---- Final solution ----
    # For methods with dense output stages (e.g., Rodas6P), use the correct solution stage.
    sol_stage = integrator.alg.num_solution_stages

    if sol_stage < stages
        # Recompute u from uprev using the solution stage row of A
        du = ks[sol_stage]
        u = uprev
        for j in 1:(sol_stage - 1)
            u = @.. u + A[sol_stage, j] * ks[j]
        end
        u = @.. u + ks[sol_stage]
    else
        # u contains uprev + Σ A[s,j]*k[j] for j=1:s-1 from the last stage loop
        du = ks[stages]
        u = @.. u + ks[stages]
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    integrator.u = u

    # ---- Error estimation ----
    if integrator.opts.adaptive
        # For Rodas methods, the error estimate is the solution stage value
        atmp = calculate_residuals(
            du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    # ---- Dense output coefficients ----
    H = tab.H
    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] = zero(integrator.k[1])
        end
        for i in 1:stages
            for j in eachindex(integrator.k)
                integrator.k[j] = @.. integrator.k[j] + H[j, i] * ks[i]
            end
        end
    end

    return nothing
end

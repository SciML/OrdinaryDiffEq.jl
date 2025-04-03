# ------------------------------------------------------------------------------
# Caches for Explicit Taylor Methods
# ------------------------------------------------------------------------------
@cache struct ExplicitTaylor2Cache{
    uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ExplicitTaylor2Cache(u, uprev, k1, k2, k3, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

struct ExplicitTaylor2ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitTaylor2ConstantCache()
end

# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylor2Cache, u) = (cache.k1, cache.k1)

@cache struct ExplicitTaylorCache{
    P, uType, rateType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    order::Val{P}
    u::uType
    uprev::uType
    us::NTuple{P, uType}
    ks::NTuple{P, rateType}
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ks = ntuple(_ -> zero(rate_prototype), Val(P))
    us = ntuple(_ -> zero(u), Val(P))
    ExplicitTaylorCache(Val(P), u, uprev, us, ks, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

struct ExplicitTaylorConstantCache{P} <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitTaylorConstantCache{P}()
end

# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylorCache, u) = (cache.ks[1], cache.ks[1])


# ------------------------------------------------------------------------------
# DAETS Caches
# ------------------------------------------------------------------------------


@cache mutable struct DAETSCache{
    uType, rateType, uNoUnitsType, tTypeNoUnits, StageLimiter, StepLimiter,
    Thread, MatType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    utilde::uType
    tmp::uType
    atmp::uType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    
    # DAE-specific
    diff_indices::Vector{Int}        # Indices for differential vs algebraic variables (0=diff, 1=alg)
    jacobian::MatType                # System Jacobian from Pantelides algorithm or something else
    taylor_coeffs::Vector{uType}     # Vector of TCs
    error_estimate::tTypeNoUnits     # Current error estimate
    u_unprojected::uType            # Unprojected solution
    residual::uType                 # Residual for projection
    correction::uType               # Correction for projection
end

function alg_cache(alg::DAETS, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, 
                   tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Type{Val{iip}}) where iip
    
    # Create standard cache components
    tmp = similar(u)
    atmp = similar(u, uEltypeNoUnits isa Type ? uEltypeNoUnits : eltype(uEltypeNoUnits))
    utilde = similar(u)
    k1 = zero(rate_prototype)
    
    # Initialize the DAE-specific components
    n = length(u)
    diff_indices = zeros(Int, n)
    for i in 1:n
        diff_indices[i] = i <= div(n, 2) ? 0 : 1
    end
    
    jacobian = Matrix{Float64}(I, n, n)
    taylor_coeffs = [similar(u) for _ in 1:10]
    u_unprojected = similar(u)
    residual = similar(u)
    correction = similar(u)
    
    # Use a concrete type for error_estimate
    error_estimate = tTypeNoUnits isa Type ? zero(tTypeNoUnits) : 0.0
    
    # Create the cache
    return DAETSCache(
        u, uprev, k1, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread,
        diff_indices, jacobian, taylor_coeffs, error_estimate,
        u_unprojected, residual, correction
    )
end

get_fsalfirstlast(cache::DAETSCache, u) = (cache.k1, cache.k1)


# ------------------------------------------------------------------------------
# DAE non-symbolic system jacobian.
# ------------------------------------------------------------------------------


"""
    build_residual_vector(dae_f, blocks, p, t)

Build the residual vector from the DAE's perspective.
In a high-index system, we might have constraints for each derivative level.
Here is a minimal example for an index-1 DAE, 
F(du, u, p, t)=0 => residual dimension = n

For a higher-order approach, you'd add additional constraints for d2u, etc.
"""
function build_residual_vector(dae_f, blocks, p, t)
    du = blocks[2]  # the first derivative block
    u  = blocks[1]
    n  = length(u)
    r  = similar(u) # residual

    # Suppose dae_f(du, u, p, t) sets the residual to 0 if consistent:
    r .= dae_f(du, u, p, t)  # dimension n

    # somehow need to add constrained for higher derivative if/when they exist

    return r
end

"""
    find_first_nonsingular_jacobian!(integrator, dae_f, u0, p, t0; max_order=5, cond_threshold=1e8)

Uses ForwardDiff to find the first non-singular Jacobian for a DAE 
    F(du, u, p, t) = 0
by successively adding higher derivatives (u', u'', etc.) if needed. 
Returns the (order, jacobian) and stores `jacobian` in the integrator's field.
"""
function find_first_nonsingular_jacobian!(
    integrator,
    dae_f::Function,
    u0::AbstractVector{<:Number},
    p,
    t0::Real;
    max_order::Int=5,
    cond_threshold::Real=1e8
)
    for order in 0:max_order
        # Build initial guess x0 = [u, 0, 0, ...]
        x0 = build_initial_guess!(u0, order)

        # Forms the residual (F=0) with up to order-th derivatives
        function res_wrapper(x::AbstractVector{<:Real})
            blocks = unpack_state_derivatives(x, u0, order)
            return build_residual_vector(dae_f, blocks, p, t0)
        end

        # Use ForwardDiff to compute Jacobian of res
        J = ForwardDiff.jacobian(res_wrapper, x0)

        if cond(J) < cond_threshold
            @info "Found well-conditioned system at order=$order"
            integrator.system_jacobian = J
            return (order, J)
        end
    end

    error("Could not find a non-singular Jacobian up to order=$max_order")
end

"""
    build_initial_guess!(u0, order)

Creates a vector [u0, 0, 0, ...] with space for state and derivatives up to order.

Returns vector of length (order+1)length(u0).
"""
function build_initial_guess!(u0, order)
    n = length(u0)
    x0 = zeros(eltype(u0), (order+1)*n)
    x0[1:n] .= u0
    return x0
end

"""
    unpack_state_derivatives(x, u0, order)

Given x = [u, du, d2u, ...], returns them as a tuple of vectors, blocks[i] being the i-th derivative.
"""
function unpack_state_derivatives(x, u0, order)
    n = length(u0)
    blocks = Vector{Vector{eltype(u0)}}(undef, order+1)
    for k in 0:order
        starti = k*n + 1
        endi   = (k+1)*n
        blocks[k+1] = @view x[starti:endi]
    end
    return blocks
end
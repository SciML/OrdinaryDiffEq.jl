using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!

# Extract the first-order derivative from TaylorScalar results.
# For scalars, access .partials[1]; for arrays, map over elements.
@inline _extract_taylor2_deriv(x::TaylorScalar) = x.partials[1]
@inline _extract_taylor2_deriv(x::AbstractArray) = map(xi -> xi.partials[1], x)

# Extract the i-th coefficient from a TaylorScalar or array of TaylorScalars.
# For scalar problems, returns a scalar. For array problems, returns a vector.
@inline _taylor_get_coefficient(ts::TaylorScalar, i::Int) = get_coefficient(ts, i)
@inline function _taylor_get_coefficient(arr::AbstractArray{<:TaylorScalar}, i::Int)
    return map(ts -> get_coefficient(ts, i), arr)
end

@inline make_taylor(all::Vararg{X, P}) where {P, X <: AbstractArray} = TaylorArray(
    Base.first(all), Base.tail(all)
)
@inline make_taylor(all::Vararg{X, P}) where {P, X} = TaylorScalar(all)

function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 2
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false
    )
    (; t, dt, uprev, u, p) = integrator
    # Unwrap FunctionWrappers since TaylorDiff types don't match wrapper signatures
    f = unwrapped_f(integrator.f)
    k1 = f(uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    k2_raw = f(u1, p, t1)
    k2 = _extract_taylor2_deriv(k2_raw)
    u = @.. uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylor2Cache)
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2Cache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; k1, k2, k3, utilde, tmp) = cache
    # Unwrap FunctionWrappers since TaylorDiff types don't match wrapper signatures
    f = unwrapped_f(integrator.f)

    # The following code is written to be fully non-allocating
    f(k1, uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    out1 = make_taylor(k1, k2)
    f(out1, u1, p, t1)
    @.. u = uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorConstantCache{P}) where {P}
    integrator.kshortsize = P
    return integrator.k = typeof(integrator.k)(undef, P)
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorConstantCache{P}, repeat_step = false
    ) where {P}
    (; t, dt, uprev, u, f, p) = integrator
    (; jet) = cache
    # jet now returns TaylorScalar with concrete numeric types
    utaylor = jet(uprev, t)
    u = eval_taylor_polynomial(utaylor, dt)
    if integrator.opts.adaptive
        utilde = TaylorDiff.get_coefficient(utaylor, P) * dt^(P + 1)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P + 1)
    # Save Taylor coefficients for dense output interpolation
    for i in 1:P
        integrator.k[i] = _taylor_get_coefficient(utaylor, i)
    end
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylorCache{P}) where {P}
    integrator.kshortsize = P
    resize!(integrator.k, P)
    # Setup k pointers
    for i in 1:P
        integrator.k[i] = get_coefficient(cache.utaylor, i)
    end
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorCache{P}, repeat_step = false
    ) where {P}
    (; t, dt, uprev, u, f, p) = integrator
    (; jet, utaylor, utilde, tmp, atmp, thread) = cache

    jet(utaylor, uprev, t)
    for i in eachindex(utaylor)
        u[i] = @inline evaluate_polynomial(utaylor[i], dt)
    end
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = TaylorDiff.get_coefficient(utaylor, P) *
            dt^(P + 1)
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P + 1)
    # Copy Taylor coefficients into k for dense output interpolation.
    # Use map! to avoid the intermediate allocation from get_coefficient.
    for i in 1:P
        map!(ts -> get_coefficient(ts, i), integrator.k[i], utaylor)
    end
    return nothing
end

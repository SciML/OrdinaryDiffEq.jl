function initialize!(integrator, cache::ExplicitRKConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::ExplicitRKConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    (; A, c, α, αEEst, stages) = cache
    (; kk) = cache

    # Calc First
    kk[1] = integrator.fsalfirst

    # Calc Middle
    for i in 2:(stages - 1)
        utilde = zero(kk[1])
        for j in 1:(i - 1)
            utilde = utilde + A[j, i] * kk[j]
        end
        kk[i] = f(uprev + dt * utilde, p, t + c[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #Calc Last
    utilde_last = zero(kk[1])
    for j in 1:(stages - 1)
        utilde_last = utilde_last + A[j, end] * kk[j]
    end
    u_beforefinal = uprev + dt * utilde_last
    kk[end] = f(u_beforefinal, p, t + c[end] * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal

    # Accumulate Result
    accum = α[1] * kk[1]
    for i in 2:stages
        accum = accum + α[i] * kk[i]
    end
    u = uprev + dt * accum

    if integrator.alg isa CompositeAlgorithm
        # Hairer II, page 22 modified to use Inf norm
        n = maximum(abs.((kk[end] .- kk[end - 1]) / (u .- u_beforefinal)))
        integrator.eigen_est = integrator.opts.internalnorm(n, t)
    end

    if integrator.opts.adaptive
        utilde = αEEst[1] .* kk[1]
        for i in 2:stages
            utilde = utilde + αEEst[i] * kk[i]
        end
        atmp = calculate_residuals(
            dt * utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if !isfsal(alg.tableau)
        integrator.fsallast = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitRKCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@generated function accumulate_explicit_stages!(
        out, A, uprev, kk, dt, ::Val{s},
        ::Val{r} = Val(s)
    ) where {s, r}
    if s == 1
        return :(@muladd @.. broadcast = false out = uprev + dt * kk[1])
    elseif s == 2
        # Note that `A` is transposed
        return :(@muladd @.. broadcast = false out = uprev + dt * (A[1, $r] * kk[1]))
    else
        expr = :(
            @muladd @.. broadcast = false out = uprev +
                dt * (A[1, $r] * kk[1] + A[2, $r] * kk[2])
        )
        acc = expr.args[end].args[end].args[end].args[end].args[end].args
        for i in 3:(s - 1)
            push!(acc, :(A[$i, $r] * kk[$i]))
        end
        return expr
    end
end

@generated function accumulate_EEst!(out, αEEst, kk, dt, ::Val{s}) where {s}
    if s == 1
        return :(@muladd @.. broadcast = false out = dt * (αEEst[1] * kk[1]))
    else
        expr = :(@muladd @.. broadcast = false out = dt * (αEEst[1] * kk[1] + αEEst[2] * kk[2]))
        acc = expr.args[end].args[end].args[end].args[end].args
        for i in 3:s
            push!(acc, :(αEEst[$i] * kk[$i]))
        end
        return expr
    end
end

function accumulate_EEst!(out, αEEst, utilde, kk, dt, stages)
    @.. broadcast = false utilde = αEEst[1] * kk[1]
    for i in 2:stages
        @.. broadcast = false utilde = utilde + αEEst[i] * kk[i]
    end
    return @.. broadcast = false out = dt * utilde
end

@muladd function compute_stages!(
        f::F, A, c, utilde, u, tmp, uprev, kk, p, t, dt,
        stages::Integer
    ) where {F}
    # Middle
    for i in 2:(stages - 1)
        @.. broadcast = false utilde = zero(kk[1][1])
        for j in 1:(i - 1)
            @.. broadcast = false utilde = utilde + A[j, i] * kk[j]
        end
        @.. broadcast = false tmp = uprev + dt * utilde
        f(kk[i], tmp, p, t + c[i] * dt)
    end

    #Last
    @.. broadcast = false utilde = zero(kk[1][1])
    for j in 1:(stages - 1)
        @.. broadcast = false utilde = utilde + A[j, end] * kk[j]
    end
    @.. broadcast = false u = uprev + dt * utilde
    f(kk[end], u, p, t + c[end] * dt) #fsallast is tmp even if not fsal
    return nothing
end

@generated function compute_stages!(
        f::F, A, c, u, tmp, uprev, kk, p, t, dt,
        ::Val{s}
    ) where {F, s}
    return quote
        Base.@nexprs $(s - 2) i′ -> begin
            i = i′ + 1
            accumulate_explicit_stages!(tmp, A, uprev, kk, dt, Val(i))
            f(kk[i], tmp, p, t + c[i] * dt)
        end
        accumulate_explicit_stages!(u, A, uprev, kk, dt, Val(s))
        f(kk[s], u, p, t + c[end] * dt)
    end
end

function runtime_split_stages!(
        f::F, A, c, utilde, u, tmp, uprev, kk, p, t, dt,
        stages::Integer
    ) where {F}
    return Base.@nif 17 (s -> (s == stages)) (
        s -> compute_stages!(
            f, A, c, u, tmp, uprev, kk, p, t,
            dt, Val(s)
        )
    ) (
        s -> compute_stages!(
            f,
            A,
            c,
            utilde,
            u,
            tmp,
            uprev,
            kk,
            p,
            t,
            dt,
            stages
        )
    )
end

function accumulate_fsal!(u, α, utilde, uprev, kk, dt, stages)
    @.. broadcast = false utilde = α[1] * kk[1]
    for i in 2:stages
        @.. broadcast = false utilde = utilde + α[i] * kk[i]
    end
    return @.. broadcast = false u = uprev + dt * utilde
end

function runtime_split_fsal!(out, A, utilde, uprev, kk, dt, stages)
    return Base.@nif 17 (s -> (s == stages)) (
        s -> accumulate_explicit_stages!(
            out, A, uprev, kk, dt,
            Val(s + 1), Val(1)
        )
    ) (
        s -> accumulate_fsal!(
            out,
            A,
            utilde,
            uprev,
            kk,
            dt,
            stages
        )
    )
end

function runtime_split_EEst!(tmp, αEEst, utilde, kk, dt, stages)
    return Base.@nif 17 (s -> (s == stages)) (s -> accumulate_EEst!(tmp, αEEst, kk, dt, Val(s))) (
        s -> accumulate_EEst!(
            tmp,
            αEEst,
            utilde,
            kk,
            dt,
            stages
        )
    )
end

@muladd function perform_step!(integrator, cache::ExplicitRKCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, false)
    # αEEst is `α - αEEst`
    (; A, c, α, αEEst, stages) = cache.tab
    (; kk, utilde, tmp, atmp) = cache

    runtime_split_stages!(f, A, c, utilde, u, tmp, uprev, kk, p, t, dt, stages)
    integrator.stats.nf += stages - 1

    #Accumulate
    if !isfsal(alg.tableau)
        runtime_split_fsal!(u, α, utilde, uprev, kk, dt, stages)
    end

    if integrator.alg isa CompositeAlgorithm
        # Hairer II, page 22 modified to use Inf norm
        @.. broadcast = false utilde = abs((kk[end] - kk[end - 1]) / (u - tmp))
        integrator.eigen_est = integrator.opts.internalnorm(norm(utilde, Inf), t)
    end

    if integrator.opts.adaptive
        runtime_split_EEst!(tmp, αEEst, utilde, kk, dt, stages)
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if !isfsal(alg.tableau)
        f(integrator.fsallast, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
end

# ============================================================================
# Generic RK interpolation utilities
# ============================================================================

@inline function falling_factorial(n::Integer, k::Integer)
    result = one(n)
    for j in 0:(k - 1)
        result *= (n - j)
    end
    return result
end

"""
    eval_poly_derivative(Θ, coeffs, order::Int)

Evaluate the k-th derivative of a polynomial at Θ using Horner's method.
coeffs[j] is the coefficient of Θ^(j-1).
"""
@inline function eval_poly_derivative(Θ, coeffs, order::Int)
    n = length(coeffs)
    if n <= order
        return zero(eltype(coeffs)) * zero(Θ)
    end

    if order == 0
        result = coeffs[n]
        for i in (n - 1):-1:1
            result = muladd(Θ, result, coeffs[i])
        end
        return result
    else
        result = coeffs[n] * falling_factorial(n - 1, order)
        for j in (n - 1):-1:(order + 1)
            multiplier = falling_factorial(j - 1, order)
            result = muladd(Θ, result, coeffs[j] * multiplier)
        end
        return result
    end
end

"""
    generic_rk_interpolant(Θ, dt, y₀, k, B_interp; idxs=nothing, order=0)

Generic interpolation for Runge-Kutta methods using the B_interp coefficient matrix.
"""
function generic_rk_interpolant(Θ, dt, y₀, k, B_interp, bi; idxs = nothing, order = 0)
    if isnothing(B_interp)
        throw(DerivativeOrderNotPossibleError())
    end

    nstages = size(B_interp, 1)

    if isempty(k) || nstages == 0
        throw(DerivativeOrderNotPossibleError())
    end

    inv_dt_factor = order <= 1 ? one(dt) : inv(dt)^(order - 1)

    # Compute weights inline (Θ may be a ForwardDiff.Dual, so we cannot
    # store into the pre-allocated Float64 buffer here).
    b1 = eval_poly_derivative(Θ, @view(B_interp[1, :]), order)
    return if isnothing(idxs)
        interp_sum = k[1] * b1
        for i in 2:nstages
            bval = eval_poly_derivative(Θ, @view(B_interp[i, :]), order)
            interp_sum = interp_sum + k[i] * bval
        end
        if order == 0
            y₀ + dt * interp_sum
        else
            interp_sum * inv_dt_factor
        end
    else
        interp_sum = k[1][idxs] * b1
        for i in 2:nstages
            bval = eval_poly_derivative(Θ, @view(B_interp[i, :]), order)
            interp_sum = interp_sum + k[i][idxs] * bval
        end
        if order == 0
            y₀[idxs] + dt * interp_sum
        else
            interp_sum * inv_dt_factor
        end
    end
end

"""
    generic_rk_interpolant!(out, Θ, dt, y₀, k, B_interp; idxs=nothing, order=0)

In-place version of generic interpolation for Runge-Kutta methods.
"""
function generic_rk_interpolant!(out, Θ, dt, y₀, k, B_interp, bi; idxs = nothing, order = 0)
    if isnothing(B_interp)
        throw(DerivativeOrderNotPossibleError())
    end

    nstages = size(B_interp, 1)

    if isempty(k) || nstages == 0
        throw(DerivativeOrderNotPossibleError())
    end

    inv_dt_factor = order <= 1 ? one(dt) : inv(dt)^(order - 1)

    # Fill pre-allocated bi buffer with polynomial weights
    for i in 1:nstages
        bi[i] = eval_poly_derivative(Θ, @view(B_interp[i, :]), order)
    end

    if isnothing(idxs)
        if order == 0
            @.. out = y₀
            for i in 1:nstages
                @.. out += dt * k[i] * bi[i]
            end
        else
            @.. out = zero(eltype(out))
            for i in 1:nstages
                @.. out += k[i] * bi[i]
            end
            @.. out *= inv_dt_factor
        end
    else
        if order == 0
            @views @.. out = y₀[idxs]
            for i in 1:nstages
                @views @.. out += dt * k[i][idxs] * bi[i]
            end
        else
            @.. out = zero(eltype(out))
            for i in 1:nstages
                @views @.. out += k[i][idxs] * bi[i]
            end
            @.. out *= inv_dt_factor
        end
    end
    return out
end

@muladd function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType, true,
        },
        integrator
    ) where {tType, uType}
    # Dispatch to the appropriate initdt algorithm
    alg_type = initdt_alg(integrator.alg)
    return _ode_determine_initdt(
        alg_type, u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob, integrator
    )
end

@muladd function _ode_determine_initdt(
        ::DefaultInitDt, u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType, true,
        },
        integrator
    ) where {tType, uType}
    _tType = eltype(tType)
    f = prob.f
    p = integrator.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    if eltype(u0) <: Number && !(integrator.alg isa CompositeAlgorithm)
        cache = get_tmp_cache(integrator)
        sk = first(cache)
        if u0 isa Array && abstol isa Number && reltol isa Number
            @inbounds @simd ivdep for i in eachindex(u0)
                sk[i] = abstol + internalnorm(u0[i], t) * reltol
            end
        else
            @.. broadcast = false sk = abstol + internalnorm(u0, t) * reltol
        end
    else
        if u0 isa Array && abstol isa Number && reltol isa Number
            sk = similar(u0, typeof(internalnorm(first(u0), t) * reltol))
            @inbounds @simd ivdep for i in eachindex(u0)
                sk[i] = abstol + internalnorm(u0[i], t) * reltol
            end
        else
            sk = @.. broadcast = false abstol + internalnorm(u0, t) * reltol
        end
    end

    if get_current_isfsal(integrator.alg, integrator.cache) &&
            integrator isa ODEIntegrator
        # Right now DelayDiffEq has issues with fsallast not being initialized
        f₀ = integrator.fsallast
        f(f₀, u0, p, t)
    else
        # TODO: use more caches
        if u0 isa Array && eltype(u0) isa Number
            T = eltype(first(u0) / oneunit_tType)
            f₀ = similar(u0, T)
            fill!(f₀, zero(T))
        else
            f₀ = zero.(u0 ./ oneunit_tType)
        end
        f(f₀, u0, p, t)
    end
    integrator.stats.nf += 1

    # TODO: use more caches
    #tmp = cache[2]

    if u0 isa Array
        if length(u0) > 0
            T = typeof(u0[1] / sk[1])
        else
            T = promote_type(eltype(u0), eltype(sk))
        end
        tmp = similar(u0, T)
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = u0[i] / sk[i]
        end
    else
        tmp = @.. broadcast = false u0 / sk
    end

    d₀ = internalnorm(tmp, t)

    #=
    Try/catch around the linear solving. This will catch singular matrices defined
    by DAEs and thus we use the _tType(1//10^(6)) default from Hairer. Note that
    this will not always catch singular matrices, an example from Andreas:

    julia> A = fill(rand(), 2, 2)
    2×2 Array{Float64,2}:
     0.637947  0.637947
     0.637947  0.637947

    julia> inv(A)
    2×2 Array{Float64,2}:
      9.0072e15  -9.0072e15
     -9.0072e15   9.0072e15

    The only way to make this more correct is to check

    issingular(A) = rank(A) < min(size(A)...)

    but that would introduce another svdfact in rank (which may not be possible
    anyways if the mass_matrix is not actually an array). Instead we stick to the
    user-chosen factorization. Sometimes this will cause `ftmp` to be absurdly
    large like shown there, but that later gets caught in the quick estimates
    below which then makes it spit out the default

    dt₀ = _tType(1//10^(6))

    so that is a a cheaper way to get to the same place than an svdfact and
    still works for matrix-free definitions of the mass matrix.
    =#

    ftmp = nothing
    if prob.f.mass_matrix != I && (
            !(prob.f isa DynamicalODEFunction) ||
                any(mm != I for mm in prob.f.mass_matrix)
        )
        ftmp = zero(f₀)
        try
            integrator.alg.linsolve(ftmp, copy(prob.f.mass_matrix), f₀, true)
            copyto!(f₀, ftmp)
        catch
            result_dt = tdir * max(smalldt, dtmin)
            @SciMLMessage(
                lazy"Mass matrix appears singular, using default small timestep: dt = $(result_dt)",
                integrator.opts.verbose, :near_singular
            )
            return result_dt
        end
    end

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = f₀[i] / sk[i] * oneunit_tType
        end
    else
        @.. broadcast = false tmp = f₀ / sk * oneunit_tType
    end

    d₁ = internalnorm(tmp, t)

    # Better than checking any(x->any(isnan, x), f₀)
    # because it also checks if partials are NaN
    # https://discourse.julialang.org/t/incorporating-forcing-functions-in-the-ode-model/70133/26
    if isnan(d₁)
        @SciMLMessage(
            "First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.",
            integrator.opts.verbose, :init_NaN
        )
        return tdir * dtmin
    end

    dt₀ = ifelse(
        (d₀ < 1 // 10^(5)) |
            (d₁ < 1 // 10^(5)), smalldt,
        convert(
            _tType,
            oneunit_tType * DiffEqBase.value(
                (d₀ / d₁) /
                    100
            )
        )
    )
    # if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    #   dt₀ = smalldt
    # else
    #   dt₀ = convert(_tType,oneunit_tType*(d₀/d₁)/100)
    # end
    dt₀ = min(dt₀, dtmax_tdir)

    if typeof(one(_tType)) <: AbstractFloat && dt₀ < 10eps(_tType) * oneunit(_tType)
        # This catches Andreas' non-singular example
        # should act like it's singular
        result_dt = tdir * max(smalldt, dtmin)
        @SciMLMessage(
            lazy"Initial timestep too small (near machine epsilon), using default: dt = $(result_dt)",
            integrator.opts.verbose, :dt_epsilon
        )
        return result_dt
    end

    dt₀_tdir = tdir * dt₀

    u₁ = zero(u0) # required by DEDataArray

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            u₁[i] = u0[i] + dt₀_tdir * f₀[i]
        end
    else
        @.. broadcast = false u₁ = u0 + dt₀_tdir * f₀
    end
    f₁ = zero(f₀)
    f(f₁, u₁, p, t + dt₀_tdir)
    integrator.stats.nf += 1

    if prob.f.mass_matrix != I && (
            !(prob.f isa DynamicalODEFunction) ||
                any(mm != I for mm in prob.f.mass_matrix)
        )
        integrator.alg.linsolve(ftmp, prob.f.mass_matrix, f₁, false)
        copyto!(f₁, ftmp)
    end

    # Constant zone before callback
    # Just return first guess
    # Avoids AD issues
    length(u0) > 0 && f₀ == f₁ && return tdir * max(dtmin, 100dt₀)

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = (f₁[i] - f₀[i]) / sk[i] * oneunit_tType
        end
    else
        @.. broadcast = false tmp = (f₁ - f₀) / sk * oneunit_tType
    end

    d₂ = internalnorm(tmp, t) / dt₀ * oneunit_tType
    # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction

    max_d₁d₂ = max(d₁, d₂)
    if max_d₁d₂ <= 1 // Int64(10)^(15)
        dt₁ = max(convert(_tType, oneunit_tType * 1 // 10^(6)), dt₀ * 1 // 10^(3))
    else
        dt₁ = convert(
            _tType,
            oneunit_tType *
                DiffEqBase.value(
                10.0^(
                    -(2 + log10(max_d₁d₂)) /
                        get_current_alg_order(
                        integrator.alg,
                        integrator.cache
                    )
                )
            )
        )
    end
    return tdir * max(dtmin, min(100dt₀, dt₁, dtmax_tdir))
end

const TYPE_NOT_CONSTANT_MESSAGE = """
Detected non-constant types in an out-of-place ODE solve, i.e. for
`du = f(u,p,t)` we see `typeof(du) !== typeof(u/t)`. This is not
supported by OrdinaryDiffEq.jl's solvers. Please either make `f`
type-constant (i.e. typeof(du) === typeof(u/t)) or use the mutating
in-place form `f(du,u,p,t)` (which is type-constant by construction).

Note that one common case for this is when computing with GPUs, using
`Float32` for `u0` and `Float64` for `tspan`. To correct this, ensure
that the element type of `tspan` matches the preferred compute type,
for example `ODEProblem(f,0f0,(0f0,1f0))` for `Float32`-based time.
"""

struct TypeNotConstantError <: Exception
    u0::Type
    f₀::Type
end

function Base.showerror(io::IO, e::TypeNotConstantError)
    println(io, TYPE_NOT_CONSTANT_MESSAGE)
    print(io, "typeof(u/t) = ")
    println(io, e.u0)
    print(io, "typeof(du) = ")
    return println(io, e.f₀)
end

@muladd function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType,
            false,
        },
        integrator
    ) where {uType, tType}
    alg_type = initdt_alg(integrator.alg)
    return _ode_determine_initdt(
        alg_type, u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob, integrator
    )
end

@muladd function _ode_determine_initdt(
        ::DefaultInitDt, u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType,
            false,
        },
        integrator
    ) where {uType, tType}
    _tType = eltype(tType)
    f = prob.f
    p = prob.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    sk = @.. broadcast = false abstol + internalnorm(u0, t) * reltol
    d₀ = internalnorm(u0 ./ sk, t)

    f₀ = f(u0, p, t)
    integrator.stats.nf += 1

    if any(x -> any(isnan, x), f₀)
        @SciMLMessage(
            "First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.",
            integrator.opts.verbose, :init_NaN
        )
    end

    inferredtype = Base.promote_op(/, typeof(u0), typeof(oneunit(t)))
    if !(f₀ isa inferredtype)
        throw(TypeNotConstantError(inferredtype, typeof(f₀)))
    end

    d₁ = internalnorm(f₀ ./ sk .* oneunit_tType, t)

    if d₀ < 1 // 10^(5) || d₁ < 1 // 10^(5)
        dt₀ = smalldt
    else
        dt₀ = convert(_tType, oneunit_tType * DiffEqBase.value((d₀ / d₁) / 100))
    end
    dt₀ = min(dt₀, dtmax_tdir)
    dt₀_tdir = tdir * dt₀

    u₁ = @.. broadcast = false u0 + dt₀_tdir * f₀
    f₁ = f(u₁, p, t + dt₀_tdir)
    integrator.stats.nf += 1

    # Constant zone before callback
    # Just return first guess
    # Avoids AD issues
    f₀ == f₁ && return tdir * max(dtmin, 100dt₀)

    d₂ = internalnorm((f₁ .- f₀) ./ sk .* oneunit_tType, t) / dt₀ * oneunit_tType

    max_d₁d₂ = max(d₁, d₂)
    if max_d₁d₂ <= 1 // Int64(10)^(15)
        dt₁ = max(smalldt, dt₀ * 1 // 10^(3))
    else
        dt₁ = _tType(
            oneunit_tType *
                DiffEqBase.value(
                10^(
                    -(2 + log10(max_d₁d₂)) /
                        get_current_alg_order(
                        integrator.alg,
                        integrator.cache
                    )
                )
            )
        )
    end
    return tdir * max(dtmin, min(100dt₀, dt₁, dtmax_tdir))
end

# Stiff initial step size algorithm (in-place)
#
# Component-wise iterative algorithm for initial step size selection, more robust
# than Hairer-Wanner for stiff multi-scale problems where some state variables
# start at zero with tiny absolute tolerances.
#
# Key differences from the Hairer-Wanner algorithm:
# 1. Uses component-wise max-norm for the upper bound (most restrictive component constrains)
# 2. Iteratively refines the estimate (up to 4 iterations)
# 3. Starts from a geometric mean of lower and upper bounds
# 4. Detects cancellation in finite-difference second derivative estimates
# 5. Applies a 0.5 safety bias factor
#
# Based on the CVHin algorithm from SUNDIALS CVODE (Hindmarsh et al., 2005).
@muladd function _ode_determine_initdt(
        ::StiffInitDt, u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType, true,
        },
        integrator
    ) where {tType, uType}
    _tType = eltype(tType)
    f = prob.f
    p = integrator.p
    oneunit_tType = oneunit(_tType)

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    tspan = prob.tspan
    tdist = abs(tspan[2] - tspan[1])

    # DAE guard: use conservative small dt for mass-matrix DAEs.
    # Must be before f₀ evaluation to avoid type issues with AD (ForwardDiff Duals).
    # Uses smalldt (≈1e-6) rather than 0.001*tdist which can be too large for
    # DAE initialization with inconsistent initial conditions.
    if integrator.isdae
        return tdir * max(smalldt, dtmin)
    end

    # Fall back to DefaultInitDt for non-Array types (GPU arrays need broadcast)
    # or non-AbstractFloat element types (ForwardDiff Duals, Complex, etc.)
    if !(u0 isa Array) || !(eltype(u0) <: AbstractFloat)
        return _ode_determine_initdt(
            DefaultInitDt(), u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob,
            integrator
        )
    end

    # Evaluate f₀ = f(u0, p, t)
    if get_current_isfsal(integrator.alg, integrator.cache) &&
            integrator isa ODEIntegrator
        f₀ = integrator.fsallast
        f(f₀, u0, p, t)
    else
        if u0 isa Array && eltype(u0) isa Number
            T = eltype(first(u0) / oneunit_tType)
            f₀ = similar(u0, T)
            fill!(f₀, zero(T))
        else
            f₀ = zero.(u0 ./ oneunit_tType)
        end
        f(f₀, u0, p, t)
    end
    integrator.stats.nf += 1

    # Handle mass matrix
    ftmp = nothing
    if prob.f.mass_matrix != I && (
            !(prob.f isa DynamicalODEFunction) ||
                any(mm != I for mm in prob.f.mass_matrix)
        )
        ftmp = zero(f₀)
        try
            integrator.alg.linsolve(ftmp, copy(prob.f.mass_matrix), f₀, true)
            copyto!(f₀, ftmp)
        catch
            result_dt = tdir * max(smalldt, dtmin)
            @SciMLMessage(
                lazy"Mass matrix appears singular, using default small timestep: dt = $(result_dt)",
                integrator.opts.verbose, :near_singular
            )
            return result_dt
        end
    end

    # Check for NaN in f₀
    if eltype(u0) <: Number && !(integrator.alg isa CompositeAlgorithm)
        cache = get_tmp_cache(integrator)
        sk = first(cache)
    else
        sk = similar(u0, typeof(internalnorm(first(u0), t) * reltol))
    end

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            sk[i] = abstol isa Number ? abstol : abstol[i]
        end
    else
        @.. broadcast = false sk = abstol
    end

    tmp = similar(u0, typeof(zero(eltype(u0)) / zero(eltype(sk)) * oneunit_tType))

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = f₀[i] / sk[i] * oneunit_tType
        end
    else
        @.. broadcast = false tmp = f₀ / sk * oneunit_tType
    end

    d₁ = internalnorm(tmp, t)

    if isnan(d₁)
        @SciMLMessage(
            "First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.",
            integrator.opts.verbose, :init_NaN
        )
        return tdir * dtmin
    end

    # === CVODE CVHin algorithm ===

    # Step 1: Compute lower and upper bounds on |h|
    # Lower bound: 100 * machine epsilon
    hlb = convert(_tType, 100 * eps(_tType) * oneunit_tType)

    # Upper bound: min of component-wise y'/y bound and 0.1 * tdist
    # hub_i = (0.1 * |y0_i| + tol_i) / |y'0_i|
    # where tol_i = rtol * |y0_i| + atol_i
    # hub = min over all i of hub_i (via max-norm of inverse)
    hub_inv = zero(_tType)
    if u0 isa Array
        @inbounds for i in eachindex(u0)
            atol_i = abstol isa Number ? abstol : abstol[i]
            rtol_i = reltol isa Number ? reltol : reltol[i]
            tol_i = rtol_i * abs(u0[i]) + atol_i
            denom = convert(_tType, 0.1) * abs(u0[i]) + tol_i
            numer = abs(f₀[i]) * oneunit_tType
            if denom > 0
                ratio = numer / denom
                hub_inv = max(hub_inv, ratio)
            end
        end
    else
        for i in eachindex(u0)
            atol_i = abstol isa Number ? abstol : abstol[i]
            rtol_i = reltol isa Number ? reltol : reltol[i]
            tol_i = rtol_i * abs(u0[i]) + atol_i
            denom = convert(_tType, 0.1) * abs(u0[i]) + tol_i
            numer = abs(f₀[i]) * oneunit_tType
            if denom > 0
                ratio = numer / denom
                hub_inv = max(hub_inv, ratio)
            end
        end
    end

    # Ensure hub_inv stays as _tType (avoid promotion from BigFloat u0 elements)
    hub_inv = convert(_tType, hub_inv)

    hub = convert(_tType, 0.1) * tdist * oneunit_tType
    if hub * hub_inv > 1
        hub = 1 / hub_inv
    end

    hub = min(hub, tdir * dtmax)

    # If bounds are crossed, return geometric mean
    if hub < hlb
        return tdir * sqrt(hlb * hub)
    end

    # Step 2: Iterative refinement
    hg = sqrt(hlb * hub)  # geometric mean = initial guess
    hs = hg               # safeguard value

    MAX_ITERS = 4
    hnew = hg

    u₁ = zero(u0)
    f₁ = zero(f₀)

    for count1 in 1:MAX_ITERS
        # Inner loop: try to evaluate ydd at trial step hg
        hg_ok = false
        for count2 in 1:MAX_ITERS
            hgs = hg * tdir

            # Euler step: u₁ = u0 + hg * f₀
            if u0 isa Array
                @inbounds @simd ivdep for i in eachindex(u0)
                    u₁[i] = u0[i] + hgs * f₀[i]
                end
            else
                @.. broadcast = false u₁ = u0 + hgs * f₀
            end

            # Evaluate f at stepped point
            f(f₁, u₁, p, t + convert(_tType, hgs))
            integrator.stats.nf += 1

            # Handle mass matrix
            if prob.f.mass_matrix != I && (
                    !(prob.f isa DynamicalODEFunction) ||
                        any(mm != I for mm in prob.f.mass_matrix)
                )
                integrator.alg.linsolve(ftmp, prob.f.mass_matrix, f₁, false)
                copyto!(f₁, ftmp)
            end

            # Check for NaN/Inf in f₁
            ydd_ok = true
            if u0 isa Array
                @inbounds for i in eachindex(f₁)
                    if !isfinite(f₁[i])
                        ydd_ok = false
                        break
                    end
                end
            else
                for i in eachindex(f₁)
                    if !isfinite(f₁[i])
                        ydd_ok = false
                        break
                    end
                end
            end

            if ydd_ok
                hg_ok = true
                break
            end

            # Reduce step and retry
            hg *= convert(_tType, 0.2)
        end

        if !hg_ok
            if count1 <= 2
                return tdir * max(smalldt, dtmin)
            end
            hnew = hs
            break
        end

        hs = hg  # save known-good step

        # Compute WRMS norm of ydd = (f₁ - f₀) / hg
        # ||ydd||_WRMS = sqrt(1/N * sum((ydd_i * ewt_i)^2))
        # where ewt_i = 1/(rtol * |y0_i| + atol_i)
        yddnrm = zero(_tType)
        N = length(u0)
        if u0 isa Array
            @inbounds for i in eachindex(u0)
                atol_i = abstol isa Number ? abstol : abstol[i]
                rtol_i = reltol isa Number ? reltol : reltol[i]
                ewt_i = 1 / (rtol_i * abs(u0[i]) + atol_i)
                ydd_i = (f₁[i] - f₀[i]) / hg * oneunit_tType
                yddnrm += (ydd_i * ewt_i)^2
            end
        else
            for i in eachindex(u0)
                atol_i = abstol isa Number ? abstol : abstol[i]
                rtol_i = reltol isa Number ? reltol : reltol[i]
                ewt_i = 1 / (rtol_i * abs(u0[i]) + atol_i)
                ydd_i = (f₁[i] - f₀[i]) / hg * oneunit_tType
                yddnrm += (ydd_i * ewt_i)^2
            end
        end
        yddnrm = convert(_tType, sqrt(yddnrm / N))

        # Compute new step proposal
        if yddnrm * hub^2 > 2
            hnew = sqrt(convert(_tType, 2) / yddnrm)
        else
            hnew = sqrt(hg * hub)
        end

        if count1 == MAX_ITERS
            break
        end

        # Convergence check
        hrat = hnew / hg
        if hrat > convert(_tType, 0.5) && hrat < 2
            break  # converged
        end

        # Cancellation detection
        if count1 > 1 && hrat > 2
            hnew = hg
            break
        end

        hg = hnew
    end

    # Step 3: Apply bias and bounds
    h0 = convert(_tType, 0.5) * hnew
    h0 = clamp(h0, hlb, hub)

    return tdir * max(dtmin, min(h0, tdir * dtmax))
end

# CVODE CVHin-style initial step size algorithm (out-of-place)
@muladd function _ode_determine_initdt(
        ::StiffInitDt, u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{
            uType, tType,
            false,
        },
        integrator
    ) where {uType, tType}
    _tType = eltype(tType)
    f = prob.f
    p = prob.p
    oneunit_tType = oneunit(_tType)

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    tspan = prob.tspan
    tdist = abs(tspan[2] - tspan[1])

    # DAE guard: use conservative small dt for mass-matrix DAEs.
    # Must be before f₀ evaluation to avoid type issues with AD (ForwardDiff Duals).
    if integrator.isdae
        return tdir * max(smalldt, dtmin)
    end

    # Fall back to DefaultInitDt for non-Array types (GPU arrays need broadcast)
    # or non-AbstractFloat element types (ForwardDiff Duals, Complex, etc.)
    if !(u0 isa Array) || !(eltype(u0) <: AbstractFloat)
        return _ode_determine_initdt(
            DefaultInitDt(), u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob,
            integrator
        )
    end

    f₀ = f(u0, p, t)
    integrator.stats.nf += 1

    if any(x -> any(isnan, x), f₀)
        @SciMLMessage(
            "First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.",
            integrator.opts.verbose, :init_NaN
        )
        return tdir * dtmin
    end

    inferredtype = Base.promote_op(/, typeof(u0), typeof(oneunit(t)))
    if !(f₀ isa inferredtype)
        throw(TypeNotConstantError(inferredtype, typeof(f₀)))
    end

    # === CVODE CVHin algorithm ===

    # Step 1: Compute lower and upper bounds on |h|
    hlb = convert(_tType, 100 * eps(_tType) * oneunit_tType)

    hub_inv = zero(_tType)
    for i in eachindex(u0)
        atol_i = abstol isa Number ? abstol : abstol[i]
        rtol_i = reltol isa Number ? reltol : reltol[i]
        tol_i = rtol_i * abs(u0[i]) + atol_i
        denom = convert(_tType, 0.1) * abs(u0[i]) + tol_i
        numer = abs(f₀[i]) * oneunit_tType
        if denom > 0
            ratio = numer / denom
            hub_inv = max(hub_inv, ratio)
        end
    end

    # Ensure hub_inv stays as _tType (avoid promotion from BigFloat u0 elements)
    hub_inv = convert(_tType, hub_inv)

    hub = convert(_tType, 0.1) * tdist * oneunit_tType
    if hub * hub_inv > 1
        hub = 1 / hub_inv
    end

    hub = min(hub, tdir * dtmax)

    if hub < hlb
        return tdir * sqrt(hlb * hub)
    end

    # Step 2: Iterative refinement
    hg = sqrt(hlb * hub)
    hs = hg

    MAX_ITERS = 4
    hnew = hg

    yddnrm = zero(_tType)

    for count1 in 1:MAX_ITERS
        hg_ok = false
        for count2 in 1:MAX_ITERS
            hgs = hg * tdir

            u₁ = @.. broadcast = false u0 + hgs * f₀
            f₁ = f(u₁, p, t + convert(_tType, hgs))
            integrator.stats.nf += 1

            ydd_ok = !any(x -> any(!isfinite, x), f₁)

            if ydd_ok
                hg_ok = true
                break
            end

            hg *= convert(_tType, 0.2)
        end

        if !hg_ok
            if count1 <= 2
                return tdir * max(smalldt, dtmin)
            end
            hnew = hs
            break
        end

        hs = hg

        hgs = hg * tdir
        u₁ = @.. broadcast = false u0 + hgs * f₀
        f₁ = f(u₁, p, t + convert(_tType, hgs))
        integrator.stats.nf += 1

        yddnrm = zero(_tType)
        N = length(u0)
        for i in eachindex(u0)
            atol_i = abstol isa Number ? abstol : abstol[i]
            rtol_i = reltol isa Number ? reltol : reltol[i]
            ewt_i = 1 / (rtol_i * abs(u0[i]) + atol_i)
            ydd_i = (f₁[i] - f₀[i]) / hg * oneunit_tType
            yddnrm += (ydd_i * ewt_i)^2
        end
        yddnrm = convert(_tType, sqrt(yddnrm / N))

        if yddnrm * hub^2 > 2
            hnew = sqrt(convert(_tType, 2) / yddnrm)
        else
            hnew = sqrt(hg * hub)
        end

        if count1 == MAX_ITERS
            break
        end

        hrat = hnew / hg
        if hrat > convert(_tType, 0.5) && hrat < 2
            break
        end

        if count1 > 1 && hrat > 2
            hnew = hg
            break
        end

        hg = hnew
    end

    # Step 3: Apply bias and bounds
    h0 = convert(_tType, 0.5) * hnew
    h0 = clamp(h0, hlb, hub)

    return tdir * max(dtmin, min(h0, tdir * dtmax))
end

# Simple initial step size for DAE problems: h = 1e-6 * tdist.
# Kept simple to avoid type issues with ComplexF64, ForwardDiff Duals, etc.
@inline function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractDAEProblem{
            duType, uType,
            tType,
        },
        integrator
    ) where {duType, uType, tType}
    _tType = eltype(tType)
    tspan = prob.tspan
    init_dt = abs(tspan[2] - tspan[1])
    init_dt = isfinite(init_dt) ? init_dt : oneunit(_tType)
    return convert(_tType, init_dt * 1 // 10^(6))
end

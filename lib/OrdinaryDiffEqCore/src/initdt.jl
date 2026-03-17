# =============================================================================
# Initial timestep estimation.
# Internal functions _ode_initdt_iip/_ode_initdt_oop implement the algorithm.
#
# When g !== nothing (stochastic), uses Hairer-Wanner with diffusion terms:
#   d₁ uses max(|f₀±3g₀|)/sk instead of f₀/sk
#   d₂ uses max(|Δf±ΔgMax|)/sk instead of Δf/sk
#
# When g === nothing (deterministic), uses CVODE CVHin algorithm
# (Hindmarsh et al., 2005) with order-dependent refinement:
#   Component-wise iterative with geometric mean of bounds,
#   step size formula h ~ (2/yddnrm)^(1/(p+1))
# =============================================================================

@muladd function _ode_initdt_iip(
        u0, t, _tType, tdir, dtmax, abstol, reltol, internalnorm,
        prob, g, noise_prototype, order, integrator
    )
    f = prob.f
    p = integrator.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    if integrator.isdae
        result_dt = tdir * max(smalldt, dtmin)
        @SciMLMessage(
            lazy"Using default small timestep for DAE: dt = $(result_dt)",
            integrator.opts.verbose, :shampine_dt
        )
        return result_dt
    end

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

    # TODO: use more caches
    #tmp = cache[2]

    if u0 isa Array
        if length(u0) > 0
            T = typeof(u0[1] / sk[1])
        else
            T = promote_type(eltype(u0), eltype(sk))
        end
        tmp = similar(u0, T)
    else
        tmp = @.. broadcast = false u0 / sk
    end

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

    if g !== nothing
        # =================================================================
        # STOCHASTIC PATH: Hairer-Wanner with diffusion terms
        # =================================================================

        if u0 isa Array
            @inbounds @simd ivdep for i in eachindex(u0)
                tmp[i] = u0[i] / sk[i]
            end
        else
            @.. broadcast = false tmp = u0 / sk
        end
        d₀ = internalnorm(tmp, t)

        # d₁: fold in diffusion terms
        # Pre-initialize g₀ so JET doesn't flag it as potentially undefined
        g₀ = nothing
        if noise_prototype !== nothing
            g₀ = zero(noise_prototype)
        else
            g₀ = zero(u0)
        end
        g(g₀, u0, p, t)
        g₀ .*= 3
        d₁ = internalnorm(
            max.(internalnorm.(f₀ .+ g₀, t), internalnorm.(f₀ .- g₀, t)) ./ sk, t
        )

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

        # d₂: fold in diffusion terms
        if noise_prototype !== nothing
            g₁ = zero(noise_prototype)
        else
            g₁ = zero(u0)
        end
        g(g₁, u₁, p, t + dt₀_tdir)
        g₁ .*= 3
        ΔgMax = max.(internalnorm.(g₀ .- g₁, t), internalnorm.(g₀ .+ g₁, t))
        d₂ = internalnorm(
            max.(
                internalnorm.(f₁ .- f₀ .+ ΔgMax, t),
                internalnorm.(f₁ .- f₀ .- ΔgMax, t)
            ) ./ sk,
            t
        ) / dt₀
        # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction

        max_d₁d₂ = max(d₁, d₂)
        if max_d₁d₂ <= 1 // Int64(10)^(15)
            dt₁ = max(convert(_tType, oneunit_tType * 1 // 10^(6)), dt₀ * 1 // 10^(3))
        else
            dt₁ = convert(
                _tType,
                oneunit_tType *
                    DiffEqBase.value(
                    10.0^(-(2 + log10(max_d₁d₂)) / order)
                )
            )
        end
        return tdir * max(dtmin, min(100dt₀, dt₁, dtmax_tdir))
    else
        # =================================================================
        # DETERMINISTIC PATH: CVHin with order dependence
        # Based on SUNDIALS CVODE's CVHin algorithm (Hindmarsh et al., 2005)
        # Component-wise iterative estimation with geometric mean of bounds
        # and order-dependent refinement: h ~ (2/yddnrm)^(1/(p+1))
        # =================================================================

        # NaN check via d₁ = norm(f₀/sk)
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

        # CVHin Step 1: Compute lower and upper bounds on |h|
        tspan = prob.tspan
        tdist = abs(tspan[2] - tspan[1])
        eps_tType = eps(_tType)
        hlb = convert(_tType, 100 * eps_tType * oneunit_tType)

        # Upper bound: most restrictive component of |f₀| / (0.1*|u0| + tol)
        hub_inv = zero(_tType)
        if u0 isa Array
            @inbounds for i in eachindex(u0)
                atol_i = abstol isa Number ? abstol : abstol[i]
                rtol_i = reltol isa Number ? reltol : reltol[i]
                tol_i = rtol_i * internalnorm(u0[i], t) + atol_i
                denom = convert(_tType, 0.1) * internalnorm(u0[i], t) + tol_i
                numer = internalnorm(f₀[i], t) * oneunit_tType
                if denom > 0
                    hub_inv = max(hub_inv, numer / denom)
                end
            end
        else
            u0_norms = internalnorm.(u0, t)
            f₀_norms = internalnorm.(f₀, t)
            tols = @.. broadcast = false reltol * u0_norms + abstol
            denoms = @.. broadcast = false convert(_tType, 0.1) * u0_norms + tols
            numers = @.. broadcast = false f₀_norms * oneunit_tType
            hub_inv_vals = ifelse.(denoms .> 0, numers ./ denoms, zero(_tType))
            hub_inv = maximum(hub_inv_vals)
        end

        hub = convert(_tType, 0.1) * tdist
        if hub * hub_inv > 1
            hub = oneunit_tType / hub_inv
        end
        hub = min(hub, abs(dtmax_tdir))

        if hub < hlb
            return tdir * max(dtmin, sqrt(hlb * hub))
        end

        # CVHin Step 2: Iterative refinement
        hg = sqrt(hlb * hub)
        hs = hg
        hnew = hg
        p_order = order
        u₁ = zero(u0)
        f₁ = zero(f₀)

        for count1 in 1:4
            # Inner loop: try hg, shrink by 0.2 if f₁ has NaN/Inf
            hg_ok = false
            for count2 in 1:4
                hgs = hg * tdir
                if u0 isa Array
                    @inbounds @simd ivdep for i in eachindex(u0)
                        u₁[i] = u0[i] + hgs * f₀[i]
                    end
                else
                    @.. broadcast = false u₁ = u0 + hgs * f₀
                end
                f(f₁, u₁, p, t + convert(_tType, hgs))

                if prob.f.mass_matrix != I && ftmp !== nothing && (
                        !(prob.f isa DynamicalODEFunction) ||
                            any(mm != I for mm in prob.f.mass_matrix)
                    )
                    integrator.alg.linsolve(ftmp, prob.f.mass_matrix, f₁, false)
                    copyto!(f₁, ftmp)
                end

                ydd_ok = true
                if u0 isa Array
                    @inbounds for i in eachindex(f₁)
                        if !isfinite(f₁[i])
                            ydd_ok = false
                            break
                        end
                    end
                else
                    ydd_ok = all(isfinite, f₁)
                end

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

            # Second derivative estimate: yddnrm = norm((f₁-f₀)/sk)/hg
            if u0 isa Array
                @inbounds @simd ivdep for i in eachindex(u0)
                    tmp[i] = (f₁[i] - f₀[i]) / sk[i] * oneunit_tType
                end
            else
                @.. broadcast = false tmp = (f₁ - f₀) / sk * oneunit_tType
            end
            yddnrm = internalnorm(tmp, t) / hg * oneunit_tType

            # Order-dependent step proposal: h ~ (2/yddnrm)^(1/(p+1))
            # Always use the formula and clamp to hub, rather than falling back
            # to the conservative geometric mean sqrt(hg*hub). This prevents
            # overly small initdt for high-order explicit methods where the
            # formula naturally gives hnew > hub.
            if DiffEqBase.value(yddnrm) > 0
                hnew = convert(
                    _tType,
                    oneunit_tType * DiffEqBase.value(
                        (2 / yddnrm)^(1 / (p_order + 1))
                    )
                )
                hnew = min(hnew, hub)
            else
                hnew = hub
            end

            count1 == 4 && break
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

        # CVHin Step 3: Apply 0.5 safety factor and bounds
        h0 = convert(_tType, 0.5) * hnew
        h0 = clamp(h0, hlb, hub)

        return tdir * max(dtmin, min(h0, abs(dtmax_tdir)))
    end
end

# ODE iip entry point
function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{uType, tType, true},
        integrator
    ) where {tType, uType}
    return _ode_initdt_iip(
        u0, t, eltype(tType), tdir, dtmax, abstol, reltol, internalnorm,
        prob, nothing, nothing,
        get_current_alg_order(integrator.alg, integrator.cache), integrator
    )
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

@muladd function _ode_initdt_oop(
        u0, t, _tType, tdir, dtmax, abstol, reltol, internalnorm,
        prob, g, order, integrator
    )
    f = prob.f
    p = prob.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(max(integrator.opts.dtmin, eps(t)))
    smalldt = max(dtmin, convert(_tType, oneunit_tType * 1 // 10^(6)))

    if integrator.isdae
        return tdir * max(smalldt, dtmin)
    end

    sk = @.. broadcast = false abstol + internalnorm(u0, t) * reltol

    f₀ = f(u0, p, t)

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

    if g !== nothing
        # =================================================================
        # STOCHASTIC PATH: Hairer-Wanner with diffusion terms
        # =================================================================

        d₀ = internalnorm(u0 ./ sk, t)

        # d₁: fold in diffusion terms
        # Pre-initialize g₀ so JET doesn't flag it as potentially undefined
        g₀ = 3g(u0, p, t)
        if any(x -> any(isnan, x), g₀)
            @SciMLMessage(
                "First function call for g produced NaNs. Exiting.",
                integrator.opts.verbose, :init_NaN
            )
        end
        d₁ = internalnorm(
            max.(internalnorm.(f₀ .+ g₀, t), internalnorm.(f₀ .- g₀, t)) ./ sk, t
        )

        if d₀ < 1 // 10^(5) || d₁ < 1 // 10^(5)
            dt₀ = smalldt
        else
            dt₀ = convert(_tType, oneunit_tType * DiffEqBase.value((d₀ / d₁) / 100))
        end
        dt₀ = min(dt₀, dtmax_tdir)
        dt₀_tdir = tdir * dt₀

        u₁ = @.. broadcast = false u0 + dt₀_tdir * f₀
        f₁ = f(u₁, p, t + dt₀_tdir)

        # Constant zone before callback
        # Just return first guess
        # Avoids AD issues
        f₀ == f₁ && return tdir * max(dtmin, 100dt₀)

        # d₂: fold in diffusion terms
        g₁ = 3g(u₁, p, t + dt₀_tdir)
        ΔgMax = max.(internalnorm.(g₀ .- g₁, t), internalnorm.(g₀ .+ g₁, t))
        d₂ = internalnorm(
            max.(
                internalnorm.(f₁ .- f₀ .+ ΔgMax, t),
                internalnorm.(f₁ .- f₀ .- ΔgMax, t)
            ) ./ sk,
            t
        ) / dt₀

        max_d₁d₂ = max(d₁, d₂)
        if max_d₁d₂ <= 1 // Int64(10)^(15)
            dt₁ = max(smalldt, dt₀ * 1 // 10^(3))
        else
            dt₁ = _tType(
                oneunit_tType *
                    DiffEqBase.value(
                    10^(-(2 + log10(max_d₁d₂)) / order)
                )
            )
        end
        return tdir * max(dtmin, min(100dt₀, dt₁, dtmax_tdir))
    else
        # =================================================================
        # DETERMINISTIC PATH: CVHin with order dependence
        # Based on SUNDIALS CVODE's CVHin algorithm (Hindmarsh et al., 2005)
        # =================================================================

        # NaN check via d₁ = norm(f₀/sk)
        d₁ = internalnorm(f₀ ./ sk .* oneunit_tType, t)
        if isnan(d₁)
            @SciMLMessage(
                "First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.",
                integrator.opts.verbose, :init_NaN
            )
            return tdir * dtmin
        end

        # CVHin Step 1: Compute lower and upper bounds
        tspan = prob.tspan
        tdist = abs(tspan[2] - tspan[1])
        eps_tType = eps(_tType)
        hlb = convert(_tType, 100 * eps_tType * oneunit_tType)

        # Upper bound: most restrictive component of |f₀| / (0.1*|u0| + tol)
        u0_norms = internalnorm.(u0, t)
        f₀_norms = internalnorm.(f₀, t)
        tols = @.. broadcast = false reltol * u0_norms + abstol
        denoms = @.. broadcast = false convert(_tType, 0.1) * u0_norms + tols
        numers = @.. broadcast = false f₀_norms * oneunit_tType
        hub_inv_vals = ifelse.(denoms .> 0, numers ./ denoms, zero(_tType))
        hub_inv = maximum(hub_inv_vals)

        hub = convert(_tType, 0.1) * tdist
        if hub * hub_inv > 1
            hub = oneunit_tType / hub_inv
        end
        hub = min(hub, abs(dtmax_tdir))

        if hub < hlb
            return tdir * max(dtmin, sqrt(hlb * hub))
        end

        # CVHin Step 2: Iterative refinement
        hg = sqrt(hlb * hub)
        hs = hg
        hnew = hg
        p_order = order

        for count1 in 1:4
            hg_ok = false
            f₁ = f₀  # will be overwritten if inner loop succeeds
            for count2 in 1:4
                hgs = hg * tdir
                u₁ = @.. broadcast = false u0 + hgs * f₀
                f₁ = f(u₁, p, t + convert(_tType, hgs))

                if !any(x -> any(!isfinite, x), f₁)
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

            # Second derivative estimate
            yddnrm = internalnorm(
                (f₁ .- f₀) ./ sk .* oneunit_tType, t
            ) / hg * oneunit_tType

            # Order-dependent step proposal: h ~ (2/yddnrm)^(1/(p+1))
            if DiffEqBase.value(yddnrm) > 0
                hnew = convert(
                    _tType,
                    oneunit_tType * DiffEqBase.value(
                        (2 / yddnrm)^(1 / (p_order + 1))
                    )
                )
                hnew = min(hnew, hub)
            else
                hnew = hub
            end

            count1 == 4 && break
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

        # CVHin Step 3: Apply 0.5 safety factor and bounds
        h0 = convert(_tType, 0.5) * hnew
        h0 = clamp(h0, hlb, hub)

        return tdir * max(dtmin, min(h0, abs(dtmax_tdir)))
    end
end

# ODE oop entry point
function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractODEProblem{uType, tType, false},
        integrator
    ) where {uType, tType}
    return _ode_initdt_oop(
        u0, t, eltype(tType), tdir, dtmax, abstol, reltol, internalnorm,
        prob, nothing,
        get_current_alg_order(integrator.alg, integrator.cache), integrator
    )
end

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

# RODE/SDE iip entry: folds noise terms into the Hairer-Wanner estimate
function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractRODEProblem{uType, tType, true},
        order, integrator
    ) where {tType, uType}
    if _get_P(integrator) !== nothing
        return tdir * dtmax / 1.0e6
    end
    g = prob.f.g
    noise_proto = hasproperty(prob, :noise_rate_prototype) ? prob.noise_rate_prototype :
        nothing
    effective_order = g !== nothing ? order + 1 // 2 : order
    return _ode_initdt_iip(
        u0, t, eltype(tType), tdir, dtmax, abstol, reltol, internalnorm,
        prob, g, noise_proto, effective_order, integrator
    )
end

# RODE/SDE oop entry: folds noise terms into the Hairer-Wanner estimate
function ode_determine_initdt(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::SciMLBase.AbstractRODEProblem{uType, tType, false},
        order, integrator
    ) where {tType, uType}
    if _get_P(integrator) !== nothing
        return tdir * dtmax / 1.0e6
    end
    g = prob.f.g
    effective_order = g !== nothing ? order + 1 // 2 : order
    return _ode_initdt_oop(
        u0, t, eltype(tType), tdir, dtmax, abstol, reltol, internalnorm,
        prob, g, effective_order, integrator
    )
end

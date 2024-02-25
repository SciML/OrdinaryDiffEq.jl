@muladd function ode_determine_initdt(u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::DiffEqBase.AbstractODEProblem{uType, tType, true
        },
        integrator) where {tType, uType}
    _tType = eltype(tType)
    f = prob.f
    p = integrator.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(integrator.opts.dtmin)
    smalldt = convert(_tType, oneunit_tType * 1 // 10^(6))

    if integrator.isdae
        return tdir * max(smalldt, dtmin)
    end

    if eltype(u0) <: Number && !(integrator.alg isa CompositeAlgorithm)
        cache = get_tmp_cache(integrator)
        sk = first(cache)
        if u0 isa Array && abstol isa Number && reltol isa Number
            @inbounds @simd ivdep for i in eachindex(u0)
                sk[i] = abstol + internalnorm(u0[i], t) * reltol
            end
        else
            @.. broadcast=false sk=abstol + internalnorm(u0, t) * reltol
        end
    else
        if u0 isa Array && abstol isa Number && reltol isa Number
            sk = similar(u0, typeof(internalnorm(first(u0), t) * reltol))
            @inbounds @simd ivdep for i in eachindex(u0)
                sk[i] = abstol + internalnorm(u0[i], t) * reltol
            end
        else
            sk = @.. broadcast=false abstol+internalnorm(u0, t) * reltol
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
            T = eltype(first(u0) / t)
            f₀ = similar(u0, T)
            fill!(f₀, zero(T))
        else
            f₀ = zero.(u0 ./ t)
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
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = u0[i] / sk[i]
        end
    else
        tmp = @.. broadcast=false u0/sk
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

    if prob.f.mass_matrix != I && (!(prob.f isa DynamicalODEFunction) ||
        any(mm != I for mm in prob.f.mass_matrix))
        ftmp = zero(f₀)
        try
            integrator.alg.linsolve(ftmp, copy(prob.f.mass_matrix), f₀, true)
            copyto!(f₀, ftmp)
        catch
            return tdir * max(smalldt, dtmin)
        end
    end

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            tmp[i] = f₀[i] / sk[i] * oneunit_tType
        end
    else
        @.. broadcast=false tmp=f₀ / sk * oneunit_tType
    end

    d₁ = internalnorm(tmp, t)

    # Better than checking any(x->any(isnan, x), f₀)
    # because it also checks if partials are NaN
    # https://discourse.julialang.org/t/incorporating-forcing-functions-in-the-ode-model/70133/26
    if integrator.opts.verbose && isnan(d₁)
        @warn("First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.")
        return tdir * dtmin
    end

    dt₀ = IfElse.ifelse((d₀ < 1 // 10^(5)) |
                        (d₁ < 1 // 10^(5)), smalldt,
        convert(_tType,
            oneunit_tType * DiffEqBase.value((d₀ / d₁) /
                                             100)))
    # if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    #   dt₀ = smalldt
    # else
    #   dt₀ = convert(_tType,oneunit_tType*(d₀/d₁)/100)
    # end
    dt₀ = min(dt₀, dtmax_tdir)

    if typeof(one(_tType)) <: AbstractFloat && dt₀ < 10eps(_tType) * oneunit(_tType)
        # This catches Andreas' non-singular example
        # should act like it's singular
        return tdir * max(smalldt, dtmin)
    end

    dt₀_tdir = tdir * dt₀

    u₁ = zero(u0) # required by DEDataArray

    if u0 isa Array
        @inbounds @simd ivdep for i in eachindex(u0)
            u₁[i] = u0[i] + dt₀_tdir * f₀[i]
        end
    else
        @.. broadcast=false u₁=u0 + dt₀_tdir * f₀
    end
    f₁ = zero(f₀)
    f(f₁, u₁, p, t + dt₀_tdir)

    if prob.f.mass_matrix != I && (!(prob.f isa DynamicalODEFunction) ||
        any(mm != I for mm in prob.f.mass_matrix))
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
        @.. broadcast=false tmp=(f₁ - f₀) / sk * oneunit_tType
    end

    d₂ = internalnorm(tmp, t) / dt₀ * oneunit_tType
    # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction

    max_d₁d₂ = max(d₁, d₂)
    if max_d₁d₂ <= 1 // Int64(10)^(15)
        dt₁ = max(convert(_tType, oneunit_tType * 1 // 10^(6)), dt₀ * 1 // 10^(3))
    else
        dt₁ = convert(_tType,
            oneunit_tType *
            DiffEqBase.value(10.0^(-(2 + log10(max_d₁d₂)) /
                                   get_current_alg_order(integrator.alg,
                integrator.cache))))
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
    println(io, e.f₀)
end

@muladd function ode_determine_initdt(u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::DiffEqBase.AbstractODEProblem{uType, tType,
            false},
        integrator) where {uType, tType}
    _tType = eltype(tType)
    f = prob.f
    p = prob.p
    oneunit_tType = oneunit(_tType)
    dtmax_tdir = tdir * dtmax

    dtmin = nextfloat(integrator.opts.dtmin)
    smalldt = convert(_tType, oneunit_tType * 1 // 10^(6))

    if integrator.isdae
        return tdir * max(smalldt, dtmin)
    end

    sk = @.. broadcast=false abstol+internalnorm(u0, t) * reltol
    d₀ = internalnorm(u0 ./ sk, t)

    f₀ = f(u0, p, t)
    if integrator.opts.verbose && any(x -> any(isnan, x), f₀)
        @warn("First function call produced NaNs. Exiting. Double check that none of the initial conditions, parameters, or timespan values are NaN.")
    end

    if Base.promote_op(/, typeof(u0), typeof(oneunit(t))) !== typeof(f₀)
        throw(TypeNotConstantError(Base.promote_op(/, typeof(u0), typeof(oneunit(t))),
            typeof(f₀)))
    end

    d₁ = internalnorm(f₀ ./ sk .* oneunit_tType, t)

    if d₀ < 1 // 10^(5) || d₁ < 1 // 10^(5)
        dt₀ = smalldt
    else
        dt₀ = convert(_tType, oneunit_tType * DiffEqBase.value((d₀ / d₁) / 100))
    end
    dt₀ = min(dt₀, dtmax_tdir)
    dt₀_tdir = tdir * dt₀

    u₁ = @.. broadcast=false u0+dt₀_tdir * f₀
    f₁ = f(u₁, p, t + dt₀_tdir)

    # Constant zone before callback
    # Just return first guess
    # Avoids AD issues
    f₀ == f₁ && return tdir * max(dtmin, 100dt₀)

    d₂ = internalnorm((f₁ .- f₀) ./ sk .* oneunit_tType, t) / dt₀ * oneunit_tType

    max_d₁d₂ = max(d₁, d₂)
    if max_d₁d₂ <= 1 // Int64(10)^(15)
        dt₁ = max(smalldt, dt₀ * 1 // 10^(3))
    else
        dt₁ = _tType(oneunit_tType *
                     DiffEqBase.value(10^(-(2 + log10(max_d₁d₂)) /
                                          get_current_alg_order(integrator.alg,
            integrator.cache))))
    end
    return tdir * max(dtmin, min(100dt₀, dt₁, dtmax_tdir))
end

@inline function ode_determine_initdt(u0, t, tdir, dtmax, abstol, reltol, internalnorm,
        prob::DiffEqBase.AbstractDAEProblem{duType, uType,
            tType},
        integrator) where {duType, uType, tType}
    _tType = eltype(tType)
    tspan = prob.tspan
    init_dt = abs(tspan[2] - tspan[1])
    init_dt = isfinite(init_dt) ? init_dt : oneunit(_tType)
    return convert(_tType, init_dt * 1 // 10^(6))
end

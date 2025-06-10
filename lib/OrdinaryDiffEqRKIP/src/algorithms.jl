using OrdinaryDiffEqCore: OrdinaryDiffEqCore

mutable struct RKIP_ALG{tableauType<:ExplicitRKTableau,elType,dtType<:AbstractVector{elType}} <: OrdinaryDiffEqAdaptiveAlgorithm
    tableau::tableauType
    dt_for_expÂ_caching::dtType
    clamp_lower_dt::Bool
    clamp_higher_dt::Bool
    use_ldiv::Bool
    cache::Union{Nothing,RKIPCache}
end

alg_order(alg::RKIP_ALG) = alg.tableau.order
alg_adaptive_order(alg::RKIP_ALG) = alg.tableau.adaptiveorder

has_dtnew_modification(alg::RKIP_ALG) = true

function dtnew_modification(alg::RKIP_ALG{tableauType,elType,dtType}, dtnew) where {tableauType,elType,dtType}
    @unpack dt_for_expÂ_caching = alg
    if first(alg.dt_for_expÂ_caching) > dtnew && alg.clamp_lower_dt
        dtnew = first(alg.dt_for_expÂ_caching)
    elseif last(alg.dt_for_expÂ_caching) < dtnew && alg.clamp_higher_dt
        dtnew = last(alg.dt_for_expÂ_caching)
    else
        dtnew = alg.dt_for_expÂ_caching[searchsortedfirst(alg.dt_for_expÂ_caching, dtnew)]
    end
    return dtnew
end

dtnew_modification(_, alg::RKIP_ALG, dtnew) = dtnew_modification(alg, dtnew)


"""
Runge-Kutta in the Interaction Picture (RIKP) method.

This is suited for solving semilinear problem of the form:

```math
\frac{du}{dt} =  Au + f(u,p,t)
```
where A is possibly stiff constant linear operator whose scaled exponential exp(Ah) can be calculated efficiently for any h.

The problem is first transformed into the Interaction Picture.

```math
u_I(t) = \exp(At) u_0
\frac{du_I}{dt} =  f_I(u_I,p,t)
f_I(u_I,p,t) = f(exp(-At)u_I, p, t)
```
This new system is then solved with an explicit Runge-Kutta method.

This type of problem is often encountered in semilinear parabolic PDE: heat diffusion, Schrödinger equation ...

This method is only implemented for semilinear problem: `SplitODEProblem` when the first function `f1` is a `AbsractSciMLOperator` A.
	
The algorithm needs for A to have the following traits of the `AbsractSciMLOperator` implemented

```julia
LinearAlgebra.exp(A, t) = #...
A(du, u, v, p, t) # for in-place problem
A(u, v p, t) # for out-of-place problem 
```

Moreover, the resulting ```expA = exp(A, t)``` must also implement ```(expA)(u, v, p, t)``` and/or ```(expA)(du, u, v, p, t)```.
The algorithm only works for constant `A` (does no depend on v, p, t).

The algorithm will cache and resue the computed operator-exponential for a set of steps. 

# Arguments
- `dtmin::T`: the smallest step `dt` for which `exp(A*dt)` will be cached. Default is `1e-3`
- `dtmax::T`: the largest step `dt` for which `exp(A*dt)` will be cached. Default is `1.0`

The set of cached steps size `dt_for_expA_caching` will follow a geometric progression such that `first(dt_for_expA_caching) = dtmin` and `last(dt_for_expA_caching) = dtmax`.

Time stepping can still happen outside the bonds but this will results in no cache (as `exp(A*h)` will be computed each step) degrading the performances. The time step can be clamped within the cache range through the args `clamp_lower_dt`  and `clamp_higher_dt`.

# Optional Arguments
- `nb_of_cache_step::Integer`: the number of steps. Default is `100`.
- `tableau::ExplicitRKTableau`: the Runge-Kutta Tableau to use, default is `constructDormandPrince6()`.
- `clamp_lower_dt::Bool`: weither to clamp proposed step to the smallest cached step in order to force the use of cached exponential, improving performance. 
	This may prevent reaching the desired tolerance. Default is `false`.
- `clamp_higher_dt::Bool`: weither to clamp proposed step to the largest cached step in order to force the use of cached exponential, improving performance. 
	This can cause performance degredation if `integrator.dtmax` is too small. Default is `true`.
- `use_ldiv::Bool`: weither, to use `ldiv(exp(A, t), v)` instead of caching `exp(A, -t)*v`. Reduce the memory usage but slightly less efficient. `ldiv` must be implemented. Only works for in-place problem. Default is `false`.

Cached operator exponentials are stored in the alorithm meaning that in the following code:

```julia
rkip = RKIP()
solve(ode_prob_1, rkip, t1)
solve(ode_prob_2, rkip, t2)
````

will reuse the cache of the first solve. This can be useful when needed to solve several time the same problem but one must be sure that `A` does not change. 

"""
RKIP(dtmin::T=1e-3, dtmax::T=1.0; nb_of_cache_step::Int=100, tableau=constructDormandPrince6(T), clamp_lower_dt::Bool=false, clamp_higher_dt::Bool=true, use_ldiv=false) where {T} = RKIP_ALG{typeof(tableau),T,Vector{T}}(
    tableau,
    logrange(dtmin, dtmax, nb_of_cache_step),
    clamp_lower_dt,
    clamp_higher_dt,
    use_ldiv,
    nothing
)

function alg_cache(alg::RKIP_ALG, u::uType, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, iip) where {uType}
    tmp = zero(u)
    utilde = zero(u)
    kk = [zero(u) for _ in 1:(alg.tableau.stages)]

    Â = isa(f, SplitFunction) ? f.f1.f : throw(ArgumentError("RKIP is only implemented for semilinear problems"))
    opType = typeof(Â)
    expOpType = typeof(exp(Â, 1.0))

    if isnothing(alg.cache)

        is_cached = Vector{Bool}(undef, length(alg.dt_for_expÂ_caching))
        fill!(is_cached, false)

        c_extended = vcat(alg.tableau.c, 1.0) # all the c values of Runge-Kutta and 1 wich is needed for the RKIP step
        c_unique = unique(c_extended) # in some tableau, there is duplicate: we only keep the unique value to save on caching time and memory
        c_index = [findfirst(==(c), c_unique) for c in c_extended] # index mapping

        exp_cache = ExpCache{expOpType}(Array{expOpType,2}(undef, length(alg.dt_for_expÂ_caching), length(c_unique)), Vector{expOpType}(undef, length(c_unique)))

        if !alg.use_ldiv
            exp_cache = ExpCacheNoLdiv(exp_cache, ExpCache{expOpType}(Array{expOpType,2}(undef, length(alg.dt_for_expÂ_caching), length(c_unique)), Vector{expOpType}(undef, length(c_unique))))
            expCacheType = ExpCacheNoLdiv{expOpType}
        else
            expCacheType = ExpCache{expOpType}
        end

        alg.cache = RKIPCache{expOpType,expCacheType,tTypeNoUnits,opType,uType,iip}(
            exp_cache,
            zero(tTypeNoUnits),
            is_cached,
            tmp,
            utilde,
            kk,
            c_unique,
            c_index
        )
    else # cache recycling 
        alg.cache = RKIPCache{expOpType,typeof(alg.cache.exp_cache),tTypeNoUnits,opType,uType,iip}(
            alg.cache.exp_cache,
            alg.last_step,
            alg.cache.is_cached,
            tmp,
            utilde,
            kk,
            alg.cache.c_unique,
            alg.cache.c_index
        )
    end
    return alg.cache
end

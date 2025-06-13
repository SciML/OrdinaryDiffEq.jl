using OrdinaryDiffEqCore: OrdinaryDiffEqCore

METHOD_DESCRIPTION = """
Runge-Kutta in the interaction picture.

This is suited for solving semilinear problem of the form:

```math
\frac{du}{dt} =  Au + f(u,p,t)
```

where A is possibly stiff time-independent linear operator whose scaled exponential exp(Ah) can be calculated efficiently for any h.
The problem is first transformed in a non-stiff variant (interaction picture)

```math
\begin{aligned}
u_I(t) &= \exp(-At) u(t) \\
\frac{du_I}{dt} &=  f_I(u_I,p,t) \\
f_I(u_I,p,t) &= f(exp(-At)u_I, p, t) \\
\end{aligned}
```
and is then solved with an explicit (adaptive) Runge-Kutta method.

This solver is only implemented for semilinear problem: `SplitODEProblem` when the first function `f1` is a `AbstractSciMLOperator` A implementing:

```julia
LinearAlgebra.exp(A, t) # = exp(A*t)
```
`A` and the return value of `exp(A, t)` must either also both implement the `AbstractSciMLOperator` interface:
```julia
A(du, u, v, p, t) # for in-place problem
A(u, v, p, t) # for out-of-place problem
```

For performance, the algorithm will cache and reuse the computed operator-exponential for a fixed set of time steps.

# Arguments
- `dtmin::T`: the smallest step `dt` for which `exp(A*dt)` will be cached. Default is `1e-3`
- `dtmax::T`: the largest step `dt` for which `exp(A*dt)` will be cached. Default is `1.0`

The fixed steps will follow a geometric progression.
Time stepping can still happen outside the bounds (for the end step for e.g) but no cache will occur (`exp(A*dt)` getting computed each step) degrading the performances.
The time step can be forcibly clamped within the cache range through the keywords `clamp_lower_dt` and `clamp_higher_dt`.

The cached operator exponentials are also directly stored in the alorithm such that:

```julia
rkip = RKIP()
solve(ode_prob_1, rkip, t1)
solve(ode_prob_2, rkip, t2)
````

will reuse the precomputed exponential cached during the first `solve` call.
This can be useful for solving several times successively problems with a common `A`.

"""
REFERENCE = """Zhongxi Zhang, Liang Chen, and Xiaoyi Bao, "A fourth-order Runge-Kutta in the interaction picture method for numerically solving the coupled nonlinear Schrödinger equation," Opt. Express 18, 8261-8276 (2010)"""

KEYWORD_DESCRIPTION = """
- `nb_of_cache_step::Integer`: the number of steps. Default is `100`.
- `tableau::ExplicitRKTableau`: the Runge-Kutta Tableau to use. Default is `constructVerner6()`.
- `clamp_lower_dt::Bool`: whether to clamp proposed step to the smallest cached step in order to force the use of cached exponential, improving performance.
	This may prevent reaching the desired tolerance. Default is `false`.
- `clamp_higher_dt::Bool`: whether to clamp proposed step to the largest cached step in order to force the use of cached exponential, improving performance.
	This can cause performance degradation if `integrator.dtmax` is too small. Default is `true`.
- `use_ldiv::Bool`: whether, to use `ldiv(exp(A, t), v)` instead of caching `exp(A, -t)*v`. Reduces the memory usage but is slightly less efficient. `ldiv` must be implemented. Only works for in-place problems. Default is `false`.
"""

@doc generic_solver_docstring(
    METHOD_DESCRIPTION, "RKIP", "Adaptative Exponential Runge-Kutta",
    REFERENCE, KEYWORD_DESCRIPTION, "")
mutable struct RKIP{
    tableauType <: ExplicitRKTableau, elType, dtType <: AbstractVector{elType}} <:
               OrdinaryDiffEqAdaptiveAlgorithm
    tableau::tableauType
    dt_for_expÂ_caching::dtType
    clamp_lower_dt::Bool
    clamp_higher_dt::Bool
    use_ldiv::Bool
    cache::Union{Nothing, RKIPCache}
end

function RKIP(dtmin::T = 1e-3, dtmax::T = 1.0; nb_of_cache_step::Int = 100,
        tableau = constructVerner6(T), clamp_lower_dt::Bool = false,
        clamp_higher_dt::Bool = true, use_ldiv = false) where {T}
    RKIP{
        typeof(tableau), T, Vector{T}}(
        tableau,
        logrange(dtmin, dtmax, nb_of_cache_step),
        clamp_lower_dt,
        clamp_higher_dt,
        use_ldiv,
        nothing
    )
end

alg_order(alg::RKIP) = alg.tableau.order
alg_adaptive_order(alg::RKIP) = alg.tableau.adaptiveorder

has_dtnew_modification(alg::RKIP) = true

function dtnew_modification(alg::RKIP{tableauType, elType, dtType},
        dtnew) where {tableauType, elType, dtType}
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

dtnew_modification(_, alg::RKIP, dtnew) = dtnew_modification(alg, dtnew)

function alg_cache(
        alg::RKIP, u::uType, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, iip) where {uType}
    tmp = zero(u)
    utilde = zero(u)
    kk = [zero(u) for _ in 1:(alg.tableau.stages)]

    Â = isa(f, SplitFunction) ? f.f1.f :
         throw(ArgumentError("RKIP is only implemented for semilinear problems"))
    opType = typeof(Â)
    expOpType = typeof(exp(Â, 1.0))

    if isnothing(alg.cache)
        is_cached = Vector{Bool}(undef, length(alg.dt_for_expÂ_caching))
        fill!(is_cached, false)

        c_extended = vcat(alg.tableau.c, 1.0) # all the c values of Runge-Kutta and 1 wich is needed for the RKIP step
        c_unique = unique(c_extended) # in some tableau, there is duplicate: we only keep the unique value to save on caching time and memory
        c_index = [findfirst(==(c), c_unique) for c in c_extended] # index mapping

        exp_cache = ExpCache{expOpType}(
            Array{expOpType, 2}(undef, length(alg.dt_for_expÂ_caching), length(c_unique)),
            Vector{expOpType}(undef, length(c_unique)))

        if !alg.use_ldiv
            exp_cache = ExpCacheNoLdiv(exp_cache,
                ExpCache{expOpType}(
                    Array{expOpType, 2}(
                        undef, length(alg.dt_for_expÂ_caching), length(c_unique)),
                    Vector{expOpType}(undef, length(c_unique))))
            expCacheType = ExpCacheNoLdiv{expOpType}
        else
            expCacheType = ExpCache{expOpType}
        end

        alg.cache = RKIPCache{expOpType, expCacheType, tTypeNoUnits, opType, uType, iip}(
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
        alg.cache = RKIPCache{
            expOpType, typeof(alg.cache.exp_cache), tTypeNoUnits, opType, uType, iip}(
            alg.cache.exp_cache,
            alg.cache.last_step,
            alg.cache.cached,
            tmp,
            utilde,
            kk,
            alg.cache.c_unique,
            alg.cache.c_mapping
        )
    end
    return alg.cache
end

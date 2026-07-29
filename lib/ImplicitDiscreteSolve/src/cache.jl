struct ImplicitDiscreteState{uType, pType, tType}
    u::uType
    p::pType
    t::tType
end

mutable struct IDSolveStepObserver{T, V}
    residualnormprev::T
    Θks::V
end

function (observer::IDSolveStepObserver)(u, fu, iteration)
    residualnorm = NonlinearSolveBase.L2_NORM(fu)
    if iteration > 1
        if observer.residualnormprev ≈ zero(observer.residualnormprev)
            push!(observer.Θks, zero(eltype(observer.Θks)))
        else
            push!(observer.Θks, residualnorm / observer.residualnormprev)
        end
    end
    observer.residualnormprev = residualnorm
    return nothing
end

function reset!(observer::IDSolveStepObserver)
    empty!(observer.Θks)
    observer.residualnormprev = zero(observer.residualnormprev)
    return nothing
end

mutable struct IDSolveCache{uType, cType, thetaType, observerType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    z::uType
    nlcache::cType
    Θks::thetaType
    observer::observerType
end

# u === nothing path: no state to evolve (e.g. MTK systems with only callbacks).
# Bypass NonlinearSolve entirely — there is nothing to nonlinear-solve, and
# constructing a NonlinearProblem with empty u would force the nonlinear
# solver to handle a degenerate case. Returning a cache with `nlcache=nothing`
# is dispatched on in `perform_step!` below to make stepping a no-op success.
function alg_cache(
        alg::IDSolve, ::Nothing, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IDSolveCache(
        nothing, nothing, nothing, nothing, uBottomEltypeNoUnits[], nothing
    )
end
function alg_cache(
        alg::IDSolve, ::Nothing, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IDSolveCache(
        nothing, nothing, nothing, nothing, uBottomEltypeNoUnits[], nothing
    )
end

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t)

    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != length(u))
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)
    Θks = uBottomEltypeNoUnits[]
    residual_prototype = isnothing(f.resid_prototype) ? u : f.resid_prototype
    observer = IDSolveStepObserver(
        zero(NonlinearSolveBase.L2_NORM(residual_prototype)), Θks
    )

    return IDSolveCache(u, uprev, state.u, nlcache, Θks, observer)
end

isdiscretecache(cache::IDSolveCache) = true

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (u_next, p) -> f(u_next, p.u, p.p, p.t)

    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != length(u))
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)
    Θks = uBottomEltypeNoUnits[]
    residual_prototype = isnothing(f.resid_prototype) ? u : f.resid_prototype
    observer = IDSolveStepObserver(
        zero(NonlinearSolveBase.L2_NORM(residual_prototype)), Θks
    )

    # FIXME Use IDSolveConstantCache?
    return IDSolveCache(u, uprev, state.u, nlcache, Θks, observer)
end

get_fsalfirstlast(cache::IDSolveCache, rate_prototype) = (nothing, nothing)

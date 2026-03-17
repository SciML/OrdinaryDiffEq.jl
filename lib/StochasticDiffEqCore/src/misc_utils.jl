struct DiffEqNLSolveTag end

struct DiffCache{T <: AbstractArray, S <: AbstractArray}
    du::T
    dual_du::S
end

Base.@pure function DiffCache(T, size, ::Type{Val{chunk_size}}) where {chunk_size}
    DiffCache(
        fill(zero(T), size...),
        fill(zero(Dual{typeof(ForwardDiff.Tag(DiffEqNLSolveTag(), T)), T, chunk_size}), size...)
    )
end

Base.@pure DiffCache(u::AbstractArray) = DiffCache(
    eltype(u), size(u), Val{ForwardDiff.pickchunksize(length(u))}
)
Base.@pure DiffCache(u::AbstractArray, nlsolve) = DiffCache(eltype(u), size(u), Val{get_chunksize(nlsolve)})
Base.@pure DiffCache(u::AbstractArray, T::Type{Val{CS}}) where {CS} = DiffCache(eltype(u), size(u), T)

get_du(dc::DiffCache, ::Type{T}) where {T <: Dual} = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

# Default nlsolve behavior, should be removed

Base.@pure determine_chunksize(u, alg::DEAlgorithm) = determine_chunksize(u, get_chunksize(alg))
Base.@pure function determine_chunksize(u, CS)
    if CS != 0
        return CS
    else
        return ForwardDiff.pickchunksize(length(u))
    end
end

struct NLSOLVEJL_SETUP{CS, AD} end
Base.@pure NLSOLVEJL_SETUP(; chunk_size = 0, autodiff = true) = NLSOLVEJL_SETUP{
    chunk_size, autodiff,
}()

# Wrapper to store the function for use with SimpleNonlinearSolve
struct IIFNLSolveFunc{F}
    f::F
end

function (p::NLSOLVEJL_SETUP{CS, AD})(f_wrapper::IIFNLSolveFunc, u0; kwargs...) where {
        CS, AD,
    }
    f = f_wrapper.f
    # Create a NonlinearProblem-compatible function
    # The IIF methods use f(resid, u) signature (in-place)
    nlf = NonlinearFunction{true}((resid, u, p) -> (f(resid, u); nothing))
    prob = NonlinearProblem(nlf, u0)
    ad = AD ? AutoForwardDiff() : AutoFiniteDiff()
    alg = SimpleTrustRegion(; autodiff = ad)
    sol = solve(prob, alg)
    return sol.u
end

function (p::NLSOLVEJL_SETUP{CS, AD})(::Type{Val{:init}}, f, u0_prototype) where {CS, AD}
    # Return a wrapper that stores the function
    return IIFNLSolveFunc(f)
end

get_chunksize(x) = 0
get_chunksize(x::NLSOLVEJL_SETUP{CS, AD}) where {CS, AD} = CS

macro cache(expr)
    name = expr.args[2].args[1].args[1]
    fields = expr.args[3].args[2:2:end]
    cache_vars = Expr[]
    rand_vars = Expr[]
    jac_vars = Pair{Symbol, Expr}[]
    ratenoise_vars = Expr[]
    for x in fields
        if x.args[2] == :uType || x.args[2] == :rateType ||
                x.args[2] == :kType || x.args[2] == :uNoUnitsType #|| x.args[2] == :possibleRateType
            push!(cache_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :JCType
            push!(cache_vars, :(c.$(x.args[1]).duals...))
        elseif x.args[2] == :GCType
            push!(cache_vars, :(c.$(x.args[1]).duals))
        elseif x.args[2] == :DiffCacheType
            push!(cache_vars, :(c.$(x.args[1]).du))
            push!(cache_vars, :(c.$(x.args[1]).dual_du))
        elseif x.args[2] == :JType || x.args[2] == :WType
            push!(jac_vars, x.args[1] => :(c.$(x.args[1])))
        elseif x.args[2] == :randType
            push!(rand_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :rateNoiseType || x.args[2] == :rateNoiseCollectionType
            # Should be a pair for handling non-diagonal
            push!(ratenoise_vars, :(c.$(x.args[1])))
        end
    end
    return quote
        $expr
        $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
        $(esc(:jac_iter))($(esc(:c))::$name) = tuple($(jac_vars...))
        $(esc(:rand_cache))($(esc(:c))::$name) = tuple($(rand_vars...))
        $(esc(:ratenoise_cache))($(esc(:c))::$name) = tuple($(ratenoise_vars...))
    end
end

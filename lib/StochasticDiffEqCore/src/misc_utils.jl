"""
    DiffEqNLSolveTag

ForwardDiff tag type used for the Jacobians of the nonlinear solves inside the SDE
integrators.

Tagging the duals keeps derivatives taken by the integrator distinguishable from
derivatives a user takes through `solve`, which is what allows nested
differentiation to work (see the ForwardDiff.jl documentation on tags).
"""
struct DiffEqNLSolveTag end

"""
    DiffCache(u)
    DiffCache(u, nlsolve)
    DiffCache(u, ::Type{Val{chunk_size}})
    DiffCache(T, size, ::Type{Val{chunk_size}})

Dual-buffered cache holding both a plain array and a `ForwardDiff.Dual` array of the
same shape.

An in-place function that must work both on ordinary numbers and on duals cannot use
a single preallocated buffer, because the element types differ. `DiffCache` stores
one buffer of each and `SciMLBase.get_du` picks the right one from the element type at
the call site.

The chunk size defaults to `ForwardDiff.pickchunksize(length(u))`, or is taken from
the nonlinear solver via [`get_chunksize`](@ref) when one is supplied.
"""
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

"""
    determine_chunksize(u, alg) -> Int
    determine_chunksize(u, CS) -> Int

ForwardDiff chunk size to use when differentiating with respect to a state like `u`.

A nonzero explicitly configured chunk size `CS` (or the one reported by
[`get_chunksize`](@ref) for `alg`) is used as-is; `0` means "not configured" and
falls back to `ForwardDiff.pickchunksize(length(u))`.
"""
Base.@pure determine_chunksize(u, alg::AbstractDEAlgorithm) = determine_chunksize(u, get_chunksize(alg))
Base.@pure function determine_chunksize(u, CS)
    if CS != 0
        return CS
    else
        return ForwardDiff.pickchunksize(length(u))
    end
end

"""
    NLSOLVEJL_SETUP(; autodiff = AutoForwardDiff())

Nonlinear-solver setup used by the IIF (implicit integrating factor) methods.

The name is historical — the solve is no longer performed by NLsolve.jl but by
`SimpleTrustRegion` from SimpleNonlinearSolve.jl. A setup object is callable in two
ways: `setup(Val{:init}, f, u0_prototype)` wraps `f` in an
[`IIFNLSolveFunc`](@ref), and `setup(wrapped_f, u0)` solves `f(resid, u) = 0`
starting from `u0` and returns the root.

## Keyword Arguments

  - `autodiff`: ADTypes.jl backend used for the Jacobian of the inner solve
    (default: `AutoForwardDiff()`).
"""
struct NLSOLVEJL_SETUP{AD}
    autodiff::AD
end
NLSOLVEJL_SETUP(; autodiff = ADTypes.AutoForwardDiff()) = NLSOLVEJL_SETUP(autodiff)

"""
    IIFNLSolveFunc(f)

Wrapper holding the residual function `f` of an IIF nonlinear solve.

[`NLSOLVEJL_SETUP`](@ref) returns one of these from its `Val{:init}` call so that the
in-place residual `f(resid, u)` can later be turned into a `NonlinearProblem` without
recapturing the surrounding cache.
"""
struct IIFNLSolveFunc{F}
    f::F
end

function (p::NLSOLVEJL_SETUP)(f_wrapper::IIFNLSolveFunc, u0; kwargs...)
    f = f_wrapper.f
    # Create a NonlinearProblem-compatible function
    # The IIF methods use f(resid, u) signature (in-place)
    nlf = NonlinearFunction{true}((resid, u, p) -> (f(resid, u); nothing))
    prob = NonlinearProblem(nlf, u0)
    alg = SimpleTrustRegion(; p.autodiff)
    sol = solve(prob, alg)
    return sol.u
end

function (p::NLSOLVEJL_SETUP)(::Type{Val{:init}}, f, u0_prototype)
    # Return a wrapper that stores the function
    return IIFNLSolveFunc(f)
end

"""
    get_chunksize(x) -> Int

ForwardDiff chunk size configured on `x`, or `0` when `x` does not configure one.

`0` is the "unset" sentinel that makes [`determine_chunksize`](@ref) fall back to
`ForwardDiff.pickchunksize`.
"""
get_chunksize(x) = 0
get_chunksize(x::NLSOLVEJL_SETUP) = OrdinaryDiffEqCore._get_fwd_chunksize_int(x.autodiff)

"""
    @cache struct MyAlgCache{...} <: StochasticDiffEqMutableCache
        ...
    end

Define an SDE solver cache and generate its buffer accessors.

`resize!` on an integrator has to grow or shrink every buffer a cache holds, and the
noise machinery has to know which of them are noise-shaped rather than state-shaped.
Rather than writing those accessors by hand for each cache, `@cache` emits them from
the declared field types:

  - `full_cache` — fields typed `uType`, `rateType`, `kType`, or `uNoUnitsType`, plus
    the `du`/`dual_du` pair of a `DiffCacheType` field and the duals of `JCType` and
    `GCType` fields.
  - `rand_cache` — fields typed `randType`.
  - `ratenoise_cache` — fields typed `rateNoiseType` or `rateNoiseCollectionType`.
  - `jac_iter` — fields typed `JType` or `WType`.

Fields whose type parameter is none of the above are left out of all four accessors,
which is the correct behavior for scalars, tableaus, and nonlinear solver objects.
"""
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
    return esc(
        quote
            $expr
            full_cache(c::$name) = tuple($(cache_vars...))
            jac_iter(c::$name) = tuple($(jac_vars...))
            rand_cache(c::$name) = tuple($(rand_vars...))
            ratenoise_cache(c::$name) = tuple($(ratenoise_vars...))
        end
    )
end

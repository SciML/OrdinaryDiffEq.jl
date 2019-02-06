abstract type AbstractNLSolver end
abstract type AbstractNLSolverCache end

mutable struct NLFunctionalCache{rateType,uType,uToltype,cType,gType} <: AbstractNLSolverCache
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  nl_iters::Int
  z::uType
  γ::gType
  c::cType
  ηold::uToltype
  # The following fields will alias for immutable cache
  z₊::uType
  dz::uType
  tmp::uType
  ztmp::uType # can be aliased with `k` if no unit
  k::rateType
end

NLFunctionalCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
  NLFunctionalCache(κ, tol, min_iter, max_iter, 0, nothing, nothing, nothing,
                    κ === nothing ? κ : zero(κ), nothing, nothing, nothing, nothing,
                    nothing)

mutable struct NLNewtonCache{rateType,uType,W,uToltype,cType,gType} <: AbstractNLSolverCache
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  nl_iters::Int
  new_W::Bool
  z::uType
  W::W
  γ::gType
  c::cType
  ηold::uToltype
  # The following fields will alias for immutable cache
  dz::uType
  tmp::uType
  ztmp::uType # can be aliased with `k` if no unit
  k::rateType
end

NLNewtonCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
  NLNewtonCache(κ, tol, min_iter, max_iter, 0, true, nothing, nothing, nothing, nothing,
                κ === nothing ? κ : zero(κ), nothing, nothing, nothing, nothing)

mutable struct NLAndersonCache{rateType,uType,uToltype,cType,gType,zsType,aType,rType} <: AbstractNLSolverCache
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  nl_iters::Int
  z::uType
  γ::gType
  c::cType
  ηold::uToltype
  alphas::aType
  residuals::rType
  # The following fields will alias for immutable cache
  z₊::uType
  dz::uType
  tmp::uType
  ztmp::uType # can be aliased with `k` if no unit
  k::rateType
  zs::zsType
  gs::zsType
end

NLAndersonCache(; κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
  NLAndersonCache(κ, tol, min_iter, max_iter, 0, nothing, nothing, nothing,
                  κ === nothing ? κ : zero(κ), ntuple(i->nothing, 9)...)

struct NLFunctional{iip,T<:NLFunctionalCache} <: AbstractNLSolver
  cache::T
end

struct NLAnderson{iip,T<:NLAndersonCache} <: AbstractNLSolver
  cache::T
  n::Int
  NLAnderson{iip,T}(nlcache::T, n=5) where {iip, T<:NLAndersonCache} = new(nlcache, n)
end

struct NLNewton{iip,T<:NLNewtonCache} <: AbstractNLSolver
  cache::T
end

# Default `iip` to `true`, but the whole type will be reinitialized in `alg_cache`
function NLFunctional(; kwargs...)
  nlcache = NLFunctionalCache(; kwargs...)
  NLFunctional{true, typeof(nlcache)}(nlcache)
end

function NLAnderson(n=5; kwargs...)
  nlcache = NLAndersonCache(; kwargs...)
  NLAnderson{true, typeof(nlcache)}(nlcache, n)
end

function NLNewton(; kwargs...)
  nlcache = NLNewtonCache(; kwargs...)
  NLNewton{true, typeof(nlcache)}(nlcache)
end

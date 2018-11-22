abstract type AbstractNLsolveSolver end
abstract type AbstractNLsolveCache end
mutable struct NLSolverCache{rateType,uType,W,uToltype,cType,gType} <: AbstractNLsolveCache
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  nl_iters::Int
  new_W::Bool
  z::uType
  W::W # NLNewton -> `W` operator; NLAnderson -> Vectors; NLFunctional -> Nothing
  γ::gType
  c::cType
  ηold::uToltype
  # The following fields will alias for immutable cache
  z₊::uType # Only used in `NLAnderson` and `NLFunctional`
  dz::uType
  tmp::uType
  b::uType # can be aliased with `k` if no unit
  k::rateType
end

struct NLFunctional{iip,T<:NLSolverCache} <: AbstractNLsolveSolver
  cache::T
end
struct NLAnderson{iip,T<:NLSolverCache} <: AbstractNLsolveSolver
  cache::T
  n::Int
end
struct NLNewton{iip,T<:NLSolverCache} <: AbstractNLsolveSolver
  cache::T
end

NLSolverCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
NLSolverCache(κ, tol, min_iter, max_iter, 0, true,
              (nothing for i in 1:10)...)

# Default `iip` to `true`, but the whole type will be reinitialized in `alg_cache`
function NLFunctional(;kwargs...)
  nlcache = NLSolverCache(;kwargs...)
  NLFunctional{true, typeof(nlcache)}(nlcache)
end
function NLAnderson(n=5; kwargs...)
  nlcache = NLSolverCache(;kwargs...)
  NLAnderson{true, typeof(nlcache)}(nlcache, n)
end
function NLNewton(;kwargs...)
  nlcache = NLSolverCache(;kwargs...)
  NLNewton{true, typeof(nlcache)}(nlcache)
end

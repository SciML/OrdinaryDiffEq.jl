abstract type AbstractNLsolveSolver end
abstract type AbstractNLsolveCache end
mutable struct NLSolverCache{rateType,uType,J,uToltype,cType,gType} <: AbstractNLsolveCache
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  new_W::Bool
  z::uType
  W::J
  γ::gType
  c::cType
  ηold::uToltype
  # The following fields will alias for mutable cache
  zprev::uType
  dz::uType
  tmp::uType
  b::uType
  k::rateType
end

struct Anderson{iip,N} <: AbstractNLsolveSolver
  cache::NLSolverCache
end
struct Newton{iip} <: AbstractNLsolveSolver
  cache::NLSolverCache
end

NLSolverCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
NLSolverCache(κ, tol, min_iter, max_iter, true,
                (nothing for i in 1:10)...)

# Default `iip` to `true`, but the whole type will be reinitialized in `alg_cache`
Anderson{N}(;kwargs...) where {N, iip} = Anderson{N,true}(NLSolverCache(kwargs...))
Functional(;kwargs...) = Anderson{1}(kwargs...)
Newton(;kwargs...) = Newton{true}(NLSolverCache(kwargs...))

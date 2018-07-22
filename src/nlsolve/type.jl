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
  # The following fields will alias for mutable cache
  zprev::uType
  dz::uType
  tmp::uType
  b::uType
  k::rateType
end

struct Functional <: AbstractNLsolveSolver
  cache::NLSolverCache
end
struct Newton <: AbstractNLsolveSolver
  cache::NLSolverCache
end

NLSolverCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
NLSolverCache(κ, tol, min_iter, max_iter, true,
              (nothing for i in 1:9)...)

Functional(;kwargs...) = Functional(NLSolverCache(kwargs...))
Newton(;kwargs...) = Functional(NLSolverCache(kwargs...))

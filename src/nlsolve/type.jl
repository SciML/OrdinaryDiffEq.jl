abstract type AbstractNLsolveSolver end
abstract type AbstractNLsolveCache end
mutable struct NLSolverCache{rateType,uType,J,uToltype,cType,gType} <: AbstractNLsolveCache
  nl_iters::Int
  κ::uToltype
  tol::uToltype
  min_iter::Int
  max_iter::Int
  new_W::Bool
  z::uType
  W::J # Newton -> `W` operator; Anderson -> Vectors; Functional -> Nothing
  γ::gType
  c::cType
  ηold::uToltype
  # The following fields will alias for immutable cache
  z₊::uType # Only used in `Anderson` and `Functional`
  dz::uType
  tmp::uType
  b::uType # can be aliased with `k` if no unit
  k::rateType
end

struct Functional{iip} <: AbstractNLsolveSolver
  cache::NLSolverCache
end
struct Anderson{iip} <: AbstractNLsolveSolver
  cache::NLSolverCache
  n::Int
end
struct Newton{iip} <: AbstractNLsolveSolver
  cache::NLSolverCache
end

NLSolverCache(;κ=nothing, tol=nothing, min_iter=1, max_iter=10) =
NLSolverCache(0, κ, tol, min_iter, max_iter, true,
              (nothing for i in 1:10)..., false)

# Default `iip` to `true`, but the whole type will be reinitialized in `alg_cache`
Functional(;kwargs...) = Functional{true}(kwargs...)
Anderson(n=5; kwargs...) = Anderson{true}(NLSolverCache(kwargs...), n)
Newton(;kwargs...) = Newton{true}(NLSolverCache(kwargs...))

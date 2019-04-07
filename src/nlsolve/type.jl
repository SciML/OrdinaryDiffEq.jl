abstract type AbstractNLSolverAlgorithm end
abstract type AbstractNLSolverCache end

@enum NLStatus::Int8 begin
  FastConvergence     = 2
  Convergence         = 1
  SlowConvergence     = 0
  VerySlowConvergence = -1
  Divergence          = -2
end

# solver

mutable struct NLSolver{iip,uType,rateType,uTolType,kType,gType,cType,C<:AbstractNLSolverCache}
  z::uType
  dz::uType
  tmp::uType
  ztmp::uType
  k::rateType
  ηold::uTolType
  κ::kType
  γ::gType
  c::cType
  max_iter::Int
  nl_iters::Int
  status::NLStatus
  cache::C
end

# algorithms

struct NLFunctional{K} <: AbstractNLSolverAlgorithm
  κ::K
  max_iter::Int
end

NLFunctional(; κ=nothing, max_iter=10) = NLFunctional(κ, max_iter)

struct NLAnderson{K,D} <: AbstractNLSolverAlgorithm
  κ::K
  max_iter::Int
  max_history::Int
  aa_start::Int
  droptol::D
end

NLAnderson(; κ=nothing, max_iter=10, max_history::Int=5, aa_start::Int=1, droptol=nothing) =
  NLAnderson(κ, max_iter, max_history, aa_start, droptol)

struct NLNewton{K} <: AbstractNLSolverAlgorithm
  κ::K
  max_iter::Int
end

NLNewton(; κ=nothing, max_iter=10) = NLNewton(κ, max_iter)

# caches

mutable struct NLNewtonCache{W,T} <: AbstractNLSolverCache
  new_W::Bool
  W::W
  Wdt::T
end

mutable struct NLNewtonConstantCache{W} <: AbstractNLSolverCache
  W::W
end

struct NLFunctionalCache{uType} <: AbstractNLSolverCache
  z₊::uType
end

struct NLFunctionalConstantCache <: AbstractNLSolverCache end

mutable struct NLAndersonCache{uType,gsType,QType,RType,gType,D} <: AbstractNLSolverCache
  z₊::uType
  dzold::uType
  z₊old::uType
  Δz₊s::gsType
  Q::QType
  R::RType
  γs::gType
  aa_start::Int
  droptol::D
end

mutable struct NLAndersonConstantCache{gsType,QType,RType,gType,D} <: AbstractNLSolverCache
  Δz₊s::gsType
  Q::QType
  R::RType
  γs::gType
  aa_start::Int
  droptol::D
end

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

mutable struct NLSolver{iip,uType,rateType,uTolType,gType,cType,C<:AbstractNLSolverCache}
  z::uType
  dz::uType
  tmp::uType
  ztmp::uType
  k::rateType
  ηold::uTolType
  κtol::uTolType
  γ::gType
  c::cType
  max_iter::Int
  nl_iters::Int
  status::NLStatus
  cache::C
end

# algorithms

struct NLFunctional{K,T} <: AbstractNLSolverAlgorithm
  κ::K
  tol::T
  max_iter::Int
end

NLFunctional(; κ=nothing, tol=nothing, max_iter=10) = NLFunctional(κ, tol, max_iter)

struct NLAnderson{K,T,D} <: AbstractNLSolverAlgorithm
  κ::K
  tol::T
  max_iter::Int
  max_history::Int
  aa_start::Int
  droptol::D
end

NLAnderson(; κ=nothing, tol=nothing, max_iter=10, max_history::Int=5, aa_start::Int=1, droptol=nothing) =
  NLAnderson(κ, tol, max_iter, max_history, aa_start, droptol)

struct NLNewton{K,T} <: AbstractNLSolverAlgorithm
  κ::K
  tol::T
  max_iter::Int
end

NLNewton(; κ=nothing, tol=nothing, max_iter=10) = NLNewton(κ, tol, max_iter)

# caches

mutable struct NLNewtonCache{W,T} <: AbstractNLSolverCache
  new_W::Bool
  W::W
  freshdt::T
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

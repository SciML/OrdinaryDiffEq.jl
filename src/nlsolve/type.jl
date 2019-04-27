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

mutable struct NLSolver{iip,uType,rateType,uTolType,kType,gType,cType,C1,C<:AbstractNLSolverCache}
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
  fast_convergence_cutoff::C1
  cache::C
end

# algorithms

struct NLFunctional{K,C} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  max_iter::Int
end

NLFunctional(; κ=nothing, max_iter=10, fast_convergence_cutoff=1//5) = NLFunctional(κ, fast_convergence_cutoff, max_iter)

struct NLAnderson{K,D,C} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  max_iter::Int
  max_history::Int
  aa_start::Int
  droptol::D
end

NLAnderson(; κ=nothing, max_iter=10, max_history::Int=5, aa_start::Int=1, droptol=nothing, fast_convergence_cutoff=1//5) =
  NLAnderson(κ, fast_convergence_cutoff, max_iter, max_history, aa_start, droptol)

struct NLNewton{K,C1,C2} <: AbstractNLSolverAlgorithm
  κ::K
  max_iter::Int
  fast_convergence_cutoff::C1
  new_W_dt_cutoff::C2
end

NLNewton(; κ=nothing, max_iter=10, fast_convergence_cutoff=1//5, new_W_dt_cutoff=1//5) = NLNewton(κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff)

# caches

mutable struct NLNewtonCache{W,T,C} <: AbstractNLSolverCache
  new_W::Bool
  W::W
  W_dt::T
  new_W_dt_cutoff::C
end

mutable struct NLNewtonConstantCache{W,C} <: AbstractNLSolverCache
  W::W
  new_W_dt_cutoff::C
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

## NLsolver

abstract type AbstractNLSolver end

mutable struct NLSolver{algType<:AbstractNLSolverAlgorithm,IIP,uType,uTolType,tTypeNoUnits,C} <: AbstractNLSolver
  """current solution"""
  z::uType
  """value `g(z)` of Newton or fixed-point iteration"""
  gz::uType
  tmp::uType
  γ::uTolType
  c::tTypeNoUnits
  alg::algType
  κ::uTolType
  η::uTolType
  fast_convergence_cutoff::uTolType
  iter::Int
  maxiters::Int
  status::NLStatus
  cache::C
end

## caches

mutable struct NLFunctionalCache{uType,tType,rateType,uNoUnitsType} <: AbstractNLSolverCache
  ustep::uType
  """residuals of `g(z) - z` of fixed-point iteration"""
  dz::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
  tstep::tType
end

mutable struct NLAndersonCache{uType,tType,rateType,uNoUnitsType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  ustep::uType
  """residuals of `g(z) - z` of fixed-point iteration"""
  dz::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  """value `g(zprev)` of previous fixed-point iteration"""
  gzprev::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzprev::uType
  Δgzs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end

mutable struct NLAndersonConstantCache{uType,tType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  """residuals `g(z) - z` of fixed-point iteration"""
  dz::uType
  tstep::tType
  """value `g(zprev)` of previous fixed-point iteration"""
  gzprev::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzprev::uType
  Δgzs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end

mutable struct NLNewtonCache{uType,tType,rateType,uNoUnitsType,J,W,du1Type,ufType,jcType,lsType,G} <: AbstractNLSolverCache
  ustep::uType
  ztmp::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  J::J
  W::W
  new_W::Bool
  W_dt::tType
  du1::du1Type
  uf::ufType
  jac_config::jcType
  linsolve::lsType
  weight::uType
  invγdt::G
  new_W_dt_cutoff::tType
end

mutable struct NLNewtonConstantCache{tType,J,W,ufType,G} <: AbstractNLSolverCache
  tstep::tType
  J::J
  W::W
  new_W::Bool
  W_dt::tType
  uf::ufType
  invγdt::G
  new_W_dt_cutoff::tType
end
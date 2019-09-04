# solver

abstract type AbstractNLSolver{algType,iip} end

mutable struct NLSolver{algType,iip,uType,uTolType,tTypeNoUnits,C<:AbstractNLSolverCache} <: AbstractNLSolver{algType,iip}
  z::uType
  tmp::uType
  ztmp::uType
  γ::uTolType
  c::tTypeNoUnits
  alg::algType
  κ::uTolType
  fast_convergence_cutoff::uTolType
  ηold::uTolType
  iter::Int
  maxiters::Int
  status::NLStatus
  cache::C
end

# caches

mutable struct NLNewtonCache{uType,tType,rateType,uNoUnitsType,J,W,du1Type,ufType,jcType,lsType,G} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  dz::uType
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
  uf::ufType
  invγdt::G
  new_W_dt_cutoff::tType
end

mutable struct NLFunctionalCache{uType,tType,rateType,uNoUnitsType} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  dz::uType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
  tstep::tType
end

mutable struct NLAndersonCache{uType,tType,rateType,uNoUnitsType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  dz::uType
  """value `g(zprev)` of previous fixed-point iteration"""
  z₊old::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzold::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::D
end

mutable struct NLAndersonConstantCache{uType,tType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  tstep::tType
  dz::uType
  """value `g(zprev)` of previous fixed-point iteration"""
  z₊old::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzold::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::D
end
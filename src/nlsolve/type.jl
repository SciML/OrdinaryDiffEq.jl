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
  iter::Int
  status::NLStatus
  fast_convergence_cutoff::C1
  cache::C
end

# caches

mutable struct NLNewtonCache{uType,tType,uNoUnitsType,J,W,du1Type,ufType,jcType,lsType,G} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
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
  uf::ufType
  invγdt::G
  new_W_dt_cutoff::tType
end

mutable struct NLFunctionalCache{uType,tType,uNoUnitsType} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  atmp::uNoUnitsType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
  tstep::tType
end

mutable struct NLAndersonCache{uType,tType,uNoUnitsType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  atmp::uNoUnitsType
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
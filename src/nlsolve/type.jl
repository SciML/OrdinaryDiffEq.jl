# solver

mutable struct NLSolver{iip,uType,rateType,uTolType,kType,gType,cType,du1Type,ufType,jcType,lsType,C1,C<:AbstractNLSolverCache}
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
  du1::du1Type
  uf::ufType
  jac_config::jcType
  linsolve::lsType
  weight::uType
  cache::C
end

# caches

mutable struct NLNewtonCache{uType,tType,uNoUnitsType,W,J,G} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  atmp::uNoUnitsType
  new_W::Bool
  W::W
  J::J
  W_dt::tType
  invγdt::G
  new_W_dt_cutoff::tType
end

mutable struct NLNewtonConstantCache{tType,W,J,G} <: AbstractNLSolverCache
  tstep::tType
  W::W
  J::J
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
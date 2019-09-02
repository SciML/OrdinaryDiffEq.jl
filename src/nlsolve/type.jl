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

mutable struct NLNewtonCache{W,J,T,C} <: AbstractNLSolverCache
  new_W::Bool
  W::W
  J::J
  W_dt::T
  new_W_dt_cutoff::C
end

mutable struct NLNewtonConstantCache{W,J,C} <: AbstractNLSolverCache
  W::W
  J::J
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

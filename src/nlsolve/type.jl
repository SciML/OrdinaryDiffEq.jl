#@enum NLStatus::Int8 begin
#  Convergence = 1
#  TryAgain = 0
#  Divergence = -1
#end

# solver

abstract type AbstractNLSolver{algType,iip} end

mutable struct NLSolver{algType,iip,uType,tType,C<:AbstractNLSolverCache} <: AbstractNLSolver{algType,iip}
  z::uType
  tmp::uType
  ztmp::uType
  γ::tType
  c::tType
  α::tType
  alg::algType
  κ::tType
  fast_convergence_cutoff::tType
  ηold::tType
  iter::Int
  maxiters::Int
  status::NLStatus
  cache::C
end

NLSolver{iip,tType}(z, tmp, ztmp, γ, c, α, alg, κ, fast_convergence_cutoff,
                    ηold, iter, maxiters, status, cache) where {iip,tType} =
  NLSolver{typeof(alg), iip, typeof(z), tType, typeof(cache)}(z, tmp, ztmp, convert(tType, γ),
                                                              convert(tType, c), convert(tType, α),
                                                              alg, convert(tType, κ),
                                                              convert(tType, fast_convergence_cutoff),
                                                              convert(tType, ηold), iter,
                                                              maxiters, status, cache)

# caches

mutable struct NLNewtonCache{uType,tType,tType2,rateType,J,W,ufType,jcType,lsType} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uType
  dz::uType
  J::J
  W::W
  new_W::Bool
  firststage::Bool
  firstcall::Bool
  W_γdt::tType
  du1::uType
  uf::ufType
  jac_config::jcType
  linsolve::lsType
  weight::uType
  invγdt::tType2
  new_W_γdt_cutoff::tType
  J_t::tType
end

mutable struct NLNewtonConstantCache{tType,tType2,J,W,ufType} <: AbstractNLSolverCache
  tstep::tType
  J::J
  W::W
  new_W::Bool
  firststage::Bool
  firstcall::Bool
  W_γdt::tType
  uf::ufType
  invγdt::tType2
  new_W_γdt_cutoff::tType
  J_t::tType
end

mutable struct NLFunctionalCache{uType,tType,rateType} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uType
  dz::uType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
  tstep::tType
end

mutable struct NLAndersonCache{uType,tType,rateType,uEltypeNoUnits} <: AbstractNLSolverCache
  ustep::uType
  tstep::tType
  k::rateType
  atmp::uType
  dz::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzold::uType
  """value `g(zprev)` of previous fixed-point iteration"""
  z₊old::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::Union{Nothing,tType}
end

mutable struct NLAndersonConstantCache{uType,tType,uEltypeNoUnits} <: AbstractNLSolverCache
  tstep::tType
  dz::uType
  """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
  dzold::uType
  """value `g(zprev)` of previous fixed-point iteration"""
  z₊old::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::Union{Nothing,tType}
end
